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
    sw = SlidingWindow(bam, min_mapq=snakemake.params["min_mapq"])  # type: ignore
    sw.write_windows(
        snakemake.output.windows,  # type: ignore
        size=750,
        step=250,
        strand_split=False,
        merge=False,
        features=True,
    )
logging.info("Done")

# # load bam file
# logging.info(f"Generating peaks from {snakemake.input['bam']}") # type: ignore
# with AlignmentFile(snakemake.input["bam"], "rb") as bam: # type: ignore
#     sw = SlidingWindow(bam, min_mapq=snakemake.params["min_mapq"]) # type: ignore
#     sw.write_windows(
#         snakemake.output.peaks, # type: ignore
#         SCHEMA,
#         size=200,
#         step=1,
#         min_rpm=2,
#         strand_split=True,
#         merge=True,
#         features=True,
#     )
# logging.info("Done")

# # close log file
# sys.stderr.close()
