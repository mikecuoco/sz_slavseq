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

params = {}
for p in snakemake.wildcards.params.split("/"):
    name, value = p.split("~")
    params[name] = int(value) if value.isdigit() else value

logger.info("Using parameters %s:", params)

# define read filter
if params["reads"] == "proper_pairs":
    read_filter = (
        lambda r: r.is_proper_pair
        and r.is_mapped
        and not r.is_supplementary
        and not r.is_secondary
        # and r.mapping_quality >= 10
    )
elif params["reads"] == "read1":
    read_filter = (
        lambda r: r.is_read1
        and r.is_mapped
        and not r.is_supplementary
        and not r.is_secondary
        # and r.mapping_quality >= 10
    )
else:
    raise Exception(f"Invalid read filter {params['reads']}")

# generate the regions
logger.info(f"Generating {params['region']} from {snakemake.input.bam}")  # type: ignore
with AlignmentFile(snakemake.input["bam"], "rb") as bam:  # type: ignore
    sw = SlidingWindow(
        bam,
        read_filter=read_filter,
        mode=params["region"],  # type: ignore
        collect_features=True if params["reads"] == "read1" else False,  # type: ignore
    )
    del params["region"]
    del params["reads"]
    try:
        sw.write_regions(
            outfile=snakemake.output[0],  # type: ignore
            min_reads=3,
            **params,
        )
        logger.info("Done")
    except:
        # import pdb; pdb.set_trace()
        # delete output file if error/interrupt
        os.remove(snakemake.output[0])  # type: ignore
        logger.error(f"Detected Error/Interuption, deleted {snakemake.output[0]}")  # type: ignore
