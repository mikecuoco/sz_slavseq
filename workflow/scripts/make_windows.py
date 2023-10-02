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
logging.info(f"Generating 750bp windows with 250bp steps from {snakemake.input['bam']}")  # type: ignore
with AlignmentFile(snakemake.input["bam"], "rb") as bam:  # type: ignore
    sw = SlidingWindow(
        bam,
        peaks=False,
        collect_features=True,
        collect_localmax=True,
    )
    sw.write_regions(
        snakemake.output.windows,  # type: ignore
        size=750,
        step=250,
        strand_split=False,
    )
logging.info("Done")

# # close log file
sys.stderr.close()
