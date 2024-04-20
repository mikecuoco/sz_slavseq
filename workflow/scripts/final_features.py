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
from pyslavseq.preprocessing import collate_labels

if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )

    logger = logging.getLogger(__name__)

    regions = pd.read_parquet(snakemake.input.regions).set_index(
        ["Chromosome", "Start", "End"]
    )

    for a in map(Path, snakemake.input.annotations):  # type: ignore
        name = a.name.split(".")[1]
        # check if file is empty
        if a.stat().st_size == 0:
            regions[name] = False
            continue
        anno = pr.read_bed(str(a)).df.set_index(["Chromosome", "Start", "End"])
        regions[name] = regions.index.isin(anno.index)
    regions.reset_index(inplace=True)

    # add en_motif scores from bigwigs
    with pyBigWig.open(snakemake.input.en_motif_pos, "r") as en_pos_bw:  # type: ignore
        regions["en_pos_score"] = regions.apply(
            lambda x: np.max(en_pos_bw.values(x["Chromosome"], x["Start"], x["End"])),
            axis=1,
        )

    with pyBigWig.open(snakemake.input.en_motif_neg, "r") as en_neg_bw:  # type: ignore
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

    regions.to_parquet(snakemake.output[0])
