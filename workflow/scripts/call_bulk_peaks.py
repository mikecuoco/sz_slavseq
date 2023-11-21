#!/usr/bin/env python
# Created on: Jul 3, 2023 at 12:21:34 AM
__author__ = "Michael Cuoco"

import sys, logging, os

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

from pysam import AlignmentFile
from pyslavseq.sliding_window import SlidingWindow
import pyarrow.parquet as pq
import pyranges as pr

logger.info(f"Generating bulk peaks from {snakemake.input.bam}")  # type: ignore
logger.info("Using parameters %s:", dict(snakemake.params))  # type: ignore
with AlignmentFile(snakemake.input.bam[0], "rb") as bam:  # type: ignore
    sw = SlidingWindow(bam, peaks=True)  # type: ignore
    try:
        sw.write_regions(
            snakemake.output.pqt,  # type: ignore
            **snakemake.params,  # type: ignore
        )
        logger.info("Done")
    except:
        # delete output file if error/interrupt
        os.remove(snakemake.output.pqt)  # type: ignore
        logger.error(f"Detected Error/Interuption, deleted {snakemake.output.pqt}")  # type: ignore
        sys.exit(1)

# convert to bed
df = pq.read_table(snakemake.output.pqt).to_pandas()  # type: ignore
pr.PyRanges(df).to_bed(snakemake.output.bed)  # type: ignore
