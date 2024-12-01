#!/usr/bin/env python
# Created on: Nov 27, 2024 at 1:19:11 PM
__author__ = "Michael Cuoco"

import numpy as np
import pandas as pd
import pyranges as pr
import pyBigWig
import logging


logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG,
)
logger = logging.getLogger(__name__)
logger.debug(f"Numpy version: {np.__version__}")
logger.debug(f"Pandas version: {pd.__version__}")
logger.debug(f"PyBigWig version: {pyBigWig.__version__}")

windows = []
with pyBigWig.open(snakemake.input.bw) as bw:
    for c in bw.chroms():
        for start, end, value in bw.intervals(c):
            windows.append(
                {"Chromosome": c, "Start": start, "End": end, "Value": value}
            )
logger.debug(f"Loaded {len(windows):,} windows")

windows.query(snakemake.params.cutoff_query, inplace=True)
logger.debug(f"Found {len(windows):,} windows passing cutoffs")

peaks = pr.PyRanges(windows).merge().df
logger.debug(f"Merged into {len(peaks):,} peaks")

# write
peaks.to_csv(snakemake.output.bed, sep="\t", index=False, header=False)
