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

# TODO: by strand?
# TODO: add 5' clipping as feature

read_filter = (
    lambda r: (not r.is_read2)
    and hasattr(r, "three_end_clipped_length")
    and (r.three_end_clipped_length > 20)
    and r.is_mapped
    and (not r.is_supplementary)
    and (not r.is_secondary)
)

# generate the regions
logger.info(f"Generating peaks from {snakemake.input.bam}")  # type: ignore
with pysam.AlignmentFile(snakemake.input.bam, "rb") as bam:  # type: ignore
    # use specific parameters for bulk samples
    gen_regions = SlidingWindow(bam, minreads=3, read_filter=read_filter).make_regions(
        collect_features=True,
        bgtest=False,
        size=5,
        step=1,
        mode="peaks",
    )

    # write to bed file
    bed_fn = snakemake.output.bed.rstrip(".gz")  # type: ignore
    cell_id = snakemake.wildcards.sample
    donor_id = snakemake.wildcards.donor
    with open(bed_fn, "w") as f:
        for i, r in enumerate(gen_regions):
            if i == 0:
                header = (
                    "#"
                    + "\t".join([str(k) for k in r.keys()])
                    + "\tcell_id\tdonor_id\n"
                )
                f.write(header)
            string = (
                "\t".join([str(v) for v in r.values()])
                + "\t"
                + cell_id
                + "\t"
                + donor_id
                + "\n"
            )
            f.write(string)

# compress and index the bed file
pysam.tabix_index(bed_fn, preset="bed", force=True)
