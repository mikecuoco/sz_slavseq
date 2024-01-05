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

from pysam import AlignmentFile
from pyslavseq.sliding_window import SlidingWindow
import pandas as pd
import pyranges as pr
import pyarrow as pa
import pyarrow.parquet as pq

params = {}
for p in snakemake.wildcards.params.split("/"):  # type: ignore
    name, value = p.split("~")
    if value.replace(".", "").isdigit():
        params[name] = int(float(value))
    elif value in ["True", "False"]:
        params[name] = True if value == "True" else False
    else:
        params[name] = value

logger.info("Using parameters %s:", params)


def write(regions: list, start):
    logger.info(
        f"Generated {len(regions)} {params['mode']} in {time.perf_counter() - start:.3f} seconds"
    )
    start = time.perf_counter()
    writer.write_table(pa.Table.from_pylist(regions))
    logger.info(
        f"Wrote {len(regions)} {params['mode']} in {time.perf_counter() - start:.3f} seconds"
    )


# generate the regions
logger.info(f"Generating {params['mode']} from {snakemake.input.bam}")  # type: ignore
with AlignmentFile(snakemake.input["bam"], "rb") as bam:  # type: ignore
    sw = SlidingWindow(bam)

    # get regions
    gen_regions = sw.make_regions(collect_features=True, **params)
    regions = [next(gen_regions)]
    schema = pa.Table.from_pylist(regions).schema
    with pq.ParquetWriter(snakemake.output.pqt, schema, write_batch_size=1e6) as writer:
        start = time.perf_counter()
        for r in gen_regions:
            if len(regions) == 1e6:
                write(regions, start)
                start = time.perf_counter()
                regions = []
            regions.append(r)

        if len(regions) > 0:
            write(regions, start)

data = pd.read_parquet(snakemake.output.pqt)
data["width"] = data["End"] - data["Start"]
pr.PyRanges(data[["Chromosome", "Start", "End", "n_reads", "width"]]).to_bed(
    snakemake.output.bed
)
