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
import numpy as np
import pyranges as pr
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyBigWig
from pyslavseq.sliding_window import SlidingWindow


# Parse parameters from parameter grid search snakemake function
def write(regions: list, start):
    """
    Write a list of regions to parquet using pyarrow
    """
    logger.info(
        f"Generated {len(regions)} peaks in {time.perf_counter() - start:.3f} seconds"
    )
    start = time.perf_counter()
    writer.write_table(pa.Table.from_pylist(regions))
    logger.info(
        f"Wrote {len(regions)} peaks in {time.perf_counter() - start:.3f} seconds"
    )


# generate the regions
logger.info(f"Generating peaks from {snakemake.input.bam}")  # type: ignore
with pysam.AlignmentFile(snakemake.input.bam, "rb") as bam:  # type: ignore
    # use specific parameters for bulk samples
    if "gDNA" in snakemake.input.bam:
        gen_regions = SlidingWindow(bam, minreads=5).make_regions(
            collect_features=True, bgtest=False, size=200, step=1, refine=True
        )
    else:
        gen_regions = SlidingWindow(bam, minreads=5).make_regions(
            collect_features=True, bgtest=True, size=200, step=1, bgsize=10000, mfold=4
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

# read back data to add final features
data = pd.read_parquet(snakemake.output.pqt)  # type: ignore
data["width"] = data["End"] - data["Start"]
warnings.filterwarnings("ignore")
pr.PyRanges(data[["Chromosome", "Start", "End", "n_reads", "width"]]).to_bed(
    snakemake.output.bed  # type: ignore
)

# # add additional features
# data["rpm"] = data["n_reads"] * (data["size_factor"] / 1e6)
# data["frac_contigs"] = data["n_contigs"] / regions["n_reads"]
# data["orientation_bias"] = (
#     np.maximum(data["n_fwd"], data["n_rev"]) / data["n_reads"]
# )
# data["frac_proper_pairs"] = data["n_proper_pairs"] / data["n_reads"]
# data["frac_duplicates"] = data["n_duplicates"] / (
#     data["n_reads"] + data["n_duplicates"]
# )
# data["frac_unique_3end"] = data["n_unique_3end"] / data["n_reads"]
# data["frac_unique_5end"] = data["n_unique_5end"] / data["n_reads"]
# data["frac_mean_supp_alignments"] = (
#     data["num_supp_alignments_mean"] / data["n_reads"]
# )

# # add en_motif scores
# with pyBigWig.open(snakemake.input.en_motif_pos, "r") as en_pos_bw  # type: ignore
#     regions["en_pos_score"] = regions.apply(
#         lambda x: np.max(en_pos_bw.values(x["Chromosome"], x["Start"], x["End"])),
#         axis=1,
#     )
# with pyBigWig.open(snakemake.input.en_motif_neg, "r") as en_neg_bw: # type: ignore
#     regions["en_neg_score"] = regions.apply(
#         lambda x: np.max(en_neg_bw.values(x["Chromosome"], x["Start"], x["End"])),
#         axis=1,
#     )

# data.to_parquet(snakemake.output.data)  # type: ignore
