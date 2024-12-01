#!/usr/bin/env python
# Created on: Oct 28, 2024 at 2:54:52 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)
from time import time
from pathlib import Path
from collections import defaultdict

import pandas as pd
from pyarrow import csv
import numpy as np
import pyranges as pr

from Bio import motifs, SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta


def compute_donor_coverage(df, bams, norm_factors):

    from pysam import AlignmentFile

    donor_id = str(
        df.name
    )  # This gets the group key (donor_id) from the groupby operation
    rff = (
        lambda x: not x.is_duplicate
        and not x.is_supplementary
        and not x.is_secondary
        and not x.is_reverse
    )
    rfr = (
        lambda x: not x.is_duplicate
        and not x.is_supplementary
        and not x.is_secondary
        and x.is_reverse
    )
    rf = lambda x: not x.is_duplicate and not x.is_supplementary and not x.is_secondary

    logger.info(f"Computing wgs coverage for {len(df)} loci for donor {donor_id}")
    # open the bams for this donor
    donor_bams = {k: AlignmentFile(v, "rb") for k, v in bams[donor_id].items()}

    def coverage(locus):
        counts = {
            "split_fwd": donor_bams["split"].count(region=locus, read_callback=rff),
            "split_rev": donor_bams["split"].count(region=locus, read_callback=rfr),
            "disc_fwd": donor_bams["disc"].count(region=locus, read_callback=rff),
            "disc_rev": donor_bams["disc"].count(region=locus, read_callback=rfr),
            "wgs": donor_bams["wgs"].count(region=locus, read_callback=rf),
        }
        return counts

    cns = df[["Chromosome", "Start", "End", "locus"]].drop_duplicates()
    cov = pd.DataFrame.from_records(
        cns["locus"].apply(coverage).tolist(), index=cns["locus"]
    )

    # close bam files
    for b in donor_bams.values():
        b.close()

    cov["split"] = cov["split_fwd"] + cov["split_rev"]
    cov["disc"] = cov["disc_fwd"] + cov["disc_rev"]
    cov["split_orientation_bias"] = (
        (cov["split_fwd"] - cov["split_rev"]).abs() / cov["split"]
    ).fillna(0)
    cov["disc_orientation_bias"] = (
        (cov["disc_fwd"] - cov["disc_rev"]).abs() / cov["disc"]
    ).fillna(0)
    for f in ["wgs", "split", "disc"]:
        cov[f"{f}_rpm"] = cov[f] / norm_factors[donor_id]

    return df.set_index("locus").join(cov, how="left").reset_index()


def label_donor(df, rmsk, mei):
    donor_id = df.name  # This gets the group key (donor_id) from the groupby operation
    genotype = load_donor_genotype(donor_id)
    anno = {**rmsk, **genotype}
    anno["MEI"] = mei
    logger.info(f"Labelling for donor {donor_id}")
    for l, d in anno.items():
        df = label(df, d, l)
    return df


if __name__ == "__main__":

    from scripts.get_peak_features import FEATURES
    from scripts.pyslavseq.io import read_tabix, load_donor_genotype, load_rmsk
    from scripts.pyslavseq.utils import label, compute_distance

    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )
    logger = logging.getLogger(__name__)

    donors = pd.read_csv(snakemake.config["donors"], sep="\t", dtype={"donor_id": str})
    libd2donor = donors.set_index("libd_id")["donor_id"].to_dict()

    # assign bams to donors
    bams = defaultdict(dict)
    for (wgs, split, disc) in zip(
        snakemake.input.wgs_bams, snakemake.input.split_bams, snakemake.input.disc_bams
    ):
        wgs, split, disc = Path(wgs), Path(split), Path(disc)
        assert (
            split.parent.name == disc.parent.name
        ), "split and disc bams must be from the same donor"
        bams[libd2donor[wgs.parent.name]]["wgs"] = wgs
        bams[split.parent.name]["split"] = split
        bams[disc.parent.name]["disc"] = disc

    norm_factors = {}
    with open(snakemake.input.wgs_coverage) as f:
        for l in f.readlines():
            libd, cov = l.split("\t")
            if libd not in libd2donor:
                continue
            norm_factors[libd2donor[libd]] = int(cov) / 1e6

    # load feature data
    data = read_tabix(snakemake.input.data)

    # compute wgs coverage
    data = (
        data.groupby("donor_id")
        .apply(
            compute_donor_coverage,
            bams=bams,
            norm_factors=norm_factors,
            include_groups=False,
        )
        .reset_index("donor_id")
    )

    # add labels
    rmsk = load_rmsk()
    megane = pr.read_bed(snakemake.input.megane).df
    graffite = pr.read_bed(snakemake.input.graffite).df
    mei = pr.PyRanges(pd.concat([megane, graffite])).extend(50).merge().df
    mei["locus"] = (
        mei["Chromosome"]
        + ":"
        + mei["Start"].astype(str)
        + "-"
        + mei["End"].astype(str)
    )
    germline = pd.concat([rmsk["L1HS"], mei])
    data = (
        data.groupby("donor_id")
        .apply(label_donor, rmsk=rmsk, mei=mei, include_groups=False)
        .reset_index("donor_id")
    )
    data = compute_distance(data, mei, "mei_distance")
    data = compute_distance(data, rmsk["L1HS"], "l1hs_distance")
    data = compute_distance(data, data, "peak_distance")
    data["germline_distance"] = data.apply(
        lambda x: min([x["mei_distance"], x["l1hs_distance"]]), axis=1
    )

    # write to pickle
    data.to_pickle(snakemake.output[0])
