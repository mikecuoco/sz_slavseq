#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys, logging
from pysam import AlignmentFile
from pyslavseq.sliding_window import SlidingWindow

# redirect stderr to log file
sys.stderr = open(snakemake.log[0], "w")  # type: ignore
logging.basicConfig(stream=sys.stderr, level=logging.INFO)

# load bam file
logging.info(f"Generating peaks from {snakemake.input['bam']}")  # type: ignore
with AlignmentFile(snakemake.input["bam"], "rb") as bam:  # type: ignore
    sw = SlidingWindow(
        bam,
        peaks=True,
        collect_features=True,
        collect_localmax=False,
    )
    sw.write_regions(
        snakemake.output.peaks,  # type: ignore
        size=150,
        step=1,
        strand_split=False,
        min_rpm=2,
        min_reads=3,
    )
logging.info("Done")

# # close log file
sys.stderr.close()
