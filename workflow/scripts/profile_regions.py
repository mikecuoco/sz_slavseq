#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import logging, os, time

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import subprocess
import pysam
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import cProfile, pstats
from pyslavseq.sliding_window import SlidingWindow

# generate the regions
logger.info(f"Profiling peak generation from {snakemake.input.bam}")  # type: ignore
with cProfile.Profile() as pr:
    with pysam.AlignmentFile(snakemake.input.bam, "rb") as bam:  # type: ignore
        gen_regions = SlidingWindow(bam).make_regions(
            collect_features=True, mode="peaks", bgtest=True, minreads=5
        )
        regions = [next(gen_regions)]
        schema = pa.Table.from_pylist(regions).schema
        with pq.ParquetWriter(
            snakemake.output.pqt, schema, write_batch_size=1e6
        ) as writer:
            for r in gen_regions:
                writer.write_table(pa.Table.from_pylist(regions))
                if len(regions) == 1e6:
                    regions = []
                regions.append(r)

            if len(regions) > 0:
                writer.write_table(pa.Table.from_pylist(regions))

with open(snakemake.output.stats, "w") as f:  # type: ignore
    p = pstats.Stats(pr, stream=f).sort_stats("cumtime")
    p.print_stats()
