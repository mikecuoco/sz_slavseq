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

# read_groups = [
#     lambda r: (r.is_read1 and r.is_forward and r.is_proper_pair) or (r.is_read2 and r.is_reverse and r.is_proper_pair),
#     lambda r: (r.is_read1 and r.is_reverse and r.is_proper_pair) or (r.is_read2 and r.is_forward and r.is_proper_pair),
# ]
read_groups = [
    lambda r: r.is_read1 and r.is_forward,
    lambda r: r.is_read1 and r.is_reverse,
]

with AlignmentFile(snakemake.input.bam[0], "rb") as bam:  # type: ignore
    for rg in read_groups:
        sw = SlidingWindow(bam, read_filter=rg, mode="peaks", collect_features=False)
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

# sort and convert to bed
df = pq.read_table(snakemake.output.pqt).to_pandas()  # type: ignore
df.sort_values(["Chromosome", "Start", "End"], inplace=True)
df.to_parquet(snakemake.output.pqt)  # type: ignore
df = df[["Chromosome", "Start", "End", "Strand"]]
pr.PyRanges(df).to_bed(snakemake.output.bed)  # type: ignore
