#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import logging, time, warnings

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import pysam
import pyranges as pr
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pyslavseq.sliding_window import SlidingWindow

params = {}
for p in snakemake.wildcards.params.split("/"):  # type: ignore
    name, value = p.split("~")
    if value.replace(".", "").isdigit():
        if "val" in name:
            params[name] = float(value)
        else:
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


# read_filter = (
#     lambda r: (not r.is_duplicate)
#     and r.is_read1
#     and r.is_proper_pair
# )

# with pysam.AlignmentFile(snakemake.input.bam, "rb") as bam:  # type: ignore
#     reads = filter(read_filter, bam.fetch())
#     tlens = [abs(r.template_length) for r in reads]

# logger.info(f"Max template length: {max(tlens)}")
# logger.info(f"Mean template length: {np.mean(tlens)}")
# logger.info(f"Median template length: {np.median(tlens)}")
# logger.info(f"Min template length: {min(tlens)}")
# size = int(np.mean(tlens) / 2)

# generate the regions
logger.info(f"Generating {params['mode']} from {snakemake.input.bam}")  # type: ignore
with pysam.AlignmentFile(snakemake.input.bam, "rb") as bam:  # type: ignore
    gen_regions = SlidingWindow(bam, minreads=5).make_regions(
        collect_features=True, **params
    )
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
warnings.filterwarnings("ignore")
pr.PyRanges(data[["Chromosome", "Start", "End", "n_reads", "width"]]).to_bed(
    snakemake.output.bed
)
