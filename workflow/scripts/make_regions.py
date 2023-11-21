#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import logging, os

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

from pysam import AlignmentFile
from pyslavseq.sliding_window import SlidingWindow

logger.info(f"Generating {snakemake.wildcards.region} from {snakemake.input.bam}")  # type: ignore
logger.info("Using parameters %s:", dict(snakemake.params))  # type: ignore

read_filter = (
    lambda r: r.is_read1
    and r.is_mapped
    and not r.is_supplementary
    and not r.is_secondary
)

with AlignmentFile(snakemake.input["bam"], "rb") as bam:  # type: ignore
    sw = SlidingWindow(
        bam,
        read_filter=read_filter,
        mode=snakemake.wildcards.region,  # type: ignore
        collect_features=True,
    )
    try:
        sw.write_regions(
            outfile=snakemake.output[0],  # type: ignore
            **snakemake.params,  # type: ignore
        )
        logger.info("Done")
    except:
        # import pdb; pdb.set_trace()
        # delete output file if error/interrupt
        os.remove(snakemake.output[0])  # type: ignore
        logger.error(f"Detected Error/Interuption, deleted {snakemake.output[0]}")  # type: ignore
