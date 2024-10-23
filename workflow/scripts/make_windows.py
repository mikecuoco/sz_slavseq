#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import logging

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import pysam
from pyslavseq.sliding_window import SlidingWindow

# generate the regions
logger.info(f"Generating windows from {snakemake.input.bam}")  # type: ignore
with pysam.AlignmentFile(snakemake.input.bam, "rb") as bam:  # type: ignore
    regions = SlidingWindow(bam, minreads=1).make_regions(
        collect_features=True,
        size=750,
        step=250,
        mode="windows",
    )

    # write to bed file
    bed_fn = snakemake.output.bed.rstrip(".gz")  # type: ignore
    with open(bed_fn, "w") as f:
        for i, r in enumerate(regions):
            if i == 0:
                header = "#" + "\t".join([str(k) for k in r.keys()]) + "\n"
                f.write(header)
            string = "\t".join([str(v) for v in r.values()]) + "\n"
            f.write(string)

# compress and index the bed file
pysam.tabix_index(bed_fn, preset="bed", force=True)
