#!/usr/bin/env python
# Created on: Jul 3, 2023 at 12:21:34 AM
__author__ = "Michael Cuoco"

import sys, logging, os

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

from pysam import AlignmentFile
from pyslavseq.sliding_window import SlidingWindow
import pyranges as pr
import pandas as pd

logger.info(f"Generating bulk peaks from {snakemake.input.bam}")  # type: ignore
logger.info("Using parameters %s:", dict(snakemake.params))  # type: ignore

with AlignmentFile(snakemake.input.bam[0], "rb") as bam:  # type: ignore
    sw = SlidingWindow(bam)
    regions = [
        r
        for r in sw.make_regions(
            mode="peaks", collect_features=True, **snakemake.params
        )
    ]

# sort and convert to bed
df = pd.DataFrame(regions).sort_values(["Chromosome", "Start", "End"])
df["width"] = df["End"] - df["Start"]
df.to_parquet(snakemake.output.pqt)  # type: ignore
df = df[["Chromosome", "Start", "End", "n_reads", "width"]]
pr.PyRanges(df).to_bed(snakemake.output.bed)  # type: ignore
