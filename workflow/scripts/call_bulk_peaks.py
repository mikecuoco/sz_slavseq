#!/usr/bin/env python
# Created on: Jul 3, 2023 at 12:21:34 AM
__author__ = "Michael Cuoco"

import sys, logging
from pysam import AlignmentFile
from pyslavseq.sliding_window import SlidingWindow
from pyslavseq.schemas import PEAKS_SCHEMA as SCHEMA

# redirect stderr to log file
sys.stderr = open(snakemake.log[0], "w")  # type: ignore
logging.basicConfig(stream=sys.stderr, level=logging.INFO)

# load bam file
with AlignmentFile(snakemake.input["bam"], "rb") as bam:  # type: ignore
    sw = SlidingWindow(bam, min_mapq=snakemake.params["min_mapq"])  # type: ignore
    sw.write_windows(
        snakemake.output.peaks,  # type: ignore
        SCHEMA,
        size=200,
        step=1,
        min_rpm=2,
        strand_split=True,
        merge=True,
        features=False,
    )

# close log file
sys.stderr.close()
