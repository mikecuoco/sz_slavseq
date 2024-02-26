#!/usr/bin/env python
# Created on: Nov 17, 2023 at 6:45:57 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

from pathlib import Path
import numpy as np
import pandas as pd
import pyranges as pr
import pyBigWig


def collate_labels(x) -> str:
    "Label one peak"
    if x["bulk_NR"] and x["megane_gaus"] and x["megane_perc"]:
        return "NONREF"
    elif x["bulk_R"] or x["ref"] or x["rmsk"] or x["primer_sites"]:
        return "REF"
    else:
        return "OTHER"


def collate_labels_bulk(x) -> str:
    "Label one peak"
    if x["bulk_NR"] or x["megane_gaus"] or x["megane_perc"]:
        return "NONREF"
    elif x["bulk_R"] or x["ref"] or x["rmsk"] or x["primer_sites"]:
        return "REF"
    else:
        return "OTHER"


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )

    logger = logging.getLogger(__name__)

    # open bigwigs
    en_pos_bw = pyBigWig.open(snakemake.input.en_motif_pos, "r")  # type: ignore
    en_neg_bw = pyBigWig.open(snakemake.input.en_motif_neg, "r")  # type: ignore

    regions = pd.read_parquet(snakemake.input.regions).set_index(
        ["Chromosome", "Start", "End"]
    )

    for a in snakemake.input.annotations:
        name = Path(a).name.split(".")[1]
        anno = pr.read_bed(a).df.set_index(["Chromosome", "Start", "End"])
        regions[name] = regions.index.isin(anno.index)

    regions.reset_index(inplace=True)

    # add en_motif scores
    regions["en_pos_score"] = regions.apply(
        lambda x: np.max(en_pos_bw.values(x["Chromosome"], x["Start"], x["End"])),
        axis=1,
    )
    regions["en_neg_score"] = regions.apply(
        lambda x: np.max(en_neg_bw.values(x["Chromosome"], x["Start"], x["End"])),
        axis=1,
    )

    # add final features
    regions["width"] = regions["End"] - regions["Start"]
    regions["rpm"] = regions["n_reads"] * (regions["size_factor"] / 1e6)
    regions["frac_contigs"] = regions["n_contigs"] / regions["n_reads"]
    regions["orientation_bias"] = (
        np.maximum(regions["n_fwd"], regions["n_rev"]) / regions["n_reads"]
    )
    regions["frac_proper_pairs"] = regions["n_proper_pairs"] / regions["n_reads"]
    regions["frac_duplicates"] = regions["n_duplicates"] / (
        regions["n_reads"] + regions["n_duplicates"]
    )
    regions["frac_unique_3end"] = regions["n_unique_3end"] / regions["n_reads"]
    regions["frac_unique_5end"] = regions["n_unique_5end"] / regions["n_reads"]
    regions["frac_mean_supp_alignments"] = (
        regions["num_supp_alignments_mean"] / regions["n_reads"]
    )

    en_pos_bw.close()
    en_neg_bw.close()

    regions.to_parquet(snakemake.output[0])
