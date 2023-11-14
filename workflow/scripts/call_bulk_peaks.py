#!/usr/bin/env python
# Created on: Jul 3, 2023 at 12:21:34 AM
__author__ = "Michael Cuoco"

import sys, logging
from pysam import AlignmentFile
from pyslavseq.sliding_window import SlidingWindow
import pyarrow.parquet as pq
import pyranges as pr

# redirect stderr to log file
sys.stderr = open(snakemake.log[0], "w")  # type: ignore
logging.basicConfig(stream=sys.stderr, level=logging.INFO)

# load bam file, save as parquet in stream
with AlignmentFile(snakemake.input.bam[0], "rb") as bam:  # type: ignore
    sw = SlidingWindow(bam, peaks=True)  # type: ignore
    sw.write_regions(
        snakemake.output.pqt,  # type: ignore
        size=200,
        step=1,
        strand_split=True,
    )

# convert to bed
df = pq.read_table(snakemake.output.pqt).to_pandas()  # type: ignore
pr.PyRanges(df).to_bed(snakemake.output.bed)  # type: ignore

# close log file
sys.stderr.close()
