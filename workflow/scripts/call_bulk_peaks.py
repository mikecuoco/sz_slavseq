#!/usr/bin/env python
# Created on: Jul 3, 2023 at 12:21:34 AM
__author__ = "Michael Cuoco"

import sys, logging
from pysam import AlignmentFile
from pyslavseq.sliding_window import SlidingWindow
from pyslavseq.schemas import PEAKS_SCHEMA as SCHEMA
import pyarrow.parquet as pq
import pyranges as pr

# redirect stderr to log file
sys.stderr = open(snakemake.log[0], "w")  # type: ignore
logging.basicConfig(stream=sys.stderr, level=logging.INFO)

# load bam file, save as parquet in stream
with AlignmentFile(snakemake.input.bam[0], "rb") as bam:  # type: ignore
    sw = SlidingWindow(bam, min_mapq=snakemake.params["min_mapq"])  # type: ignore
    sw.write_windows(
        snakemake.output.pqt,  # type: ignore
        SCHEMA,
        size=200,
        step=1,
        strand_split=True,
        merge=True,
        features=False,
    )

# convert to bed
df = pq.read_table(snakemake.output.pqt).to_pandas()  # type: ignore
pr.PyRanges(df).to_bed(snakemake.output.bed)  # type: ignore

# close log file
sys.stderr.close()
