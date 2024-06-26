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
    if "gDNA" in snakemake.input.bam:  # type: ignore
        gen_regions = SlidingWindow(bam, minreads=5).make_regions(
            collect_features=True, bgtest=False, size=200, step=1
        )
    else:
        gen_regions = SlidingWindow(bam, minreads=5).make_regions(
            collect_features=True, bgtest=True, size=200, step=1, bgsize=10000, mfold=4
        )
    regions = [next(gen_regions)]
    schema = pa.Table.from_pylist(regions).schema
    with pq.ParquetWriter(snakemake.output.pqt, schema, write_batch_size=1e6) as writer:  # type: ignore
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

# validate
for i, r in data.iterrows():
    for bg in [5e3, 1e4, 2e4]:
        assert (
            r[f"bg{int(bg)}_Chromosome"] == r["Chromosome"]
        ), "Some chromosomes do not match"
        assert r[f"bg{int(bg)}_Start"] <= r["Start"], "Some starts do not match"
        assert r[f"bg{int(bg)}_End"] >= r["End"], "Some ends do not match"

# # add additional features
for p in ["", "bg5000_", "bg10000_", "bg20000_"]:
    data[f"{p}rpm"] = data[f"{p}n_reads"] * (data["size_factor"] / 1e6)
    data[f"{p}frac_contigs"] = data[f"{p}n_contigs"] / data[f"{p}n_reads"]
    data[f"{p}orientation_bias"] = (
        np.maximum(data[f"{p}n_fwd"], data[f"{p}n_rev"]) / data[f"{p}n_reads"]
    )
    data[f"{p}frac_proper_pairs"] = data[f"{p}n_proper_pairs"] / data[f"{p}n_reads"]
    data[f"{p}frac_duplicates"] = data[f"{p}n_duplicates"] / (
        data[f"{p}n_reads"] + data[f"{p}n_duplicates"]
    )
    data[f"{p}frac_unique_3end"] = data[f"{p}n_unique_3end"] / data[f"{p}n_reads"]
    data[f"{p}frac_unique_5end"] = data[f"{p}n_unique_5end"] / data[f"{p}n_reads"]
    data[f"{p}frac_mean_supp_alignments"] = (
        data[f"{p}num_supp_alignments_mean"] / data[f"{p}n_reads"]
    )

    # # add en_motif scores
    with pyBigWig.open(snakemake.input.en_motif_pos, "r") as en_pos_bw:  # type: ignore
        data[f"{p}en_pos_score"] = data.apply(
            lambda x: np.max(
                en_pos_bw.values(x[f"{p}Chromosome"], x[f"{p}Start"], x[f"{p}End"])
            ),
            axis=1,
        )
    with pyBigWig.open(snakemake.input.en_motif_neg, "r") as en_neg_bw:  # type: ignore
        data[f"{p}en_neg_score"] = data.apply(
            lambda x: np.max(
                en_neg_bw.values(x[f"{p}Chromosome"], x[f"{p}Start"], x[f"{p}End"])
            ),
            axis=1,
        )

data.to_parquet(snakemake.output.pqt)  # type: ignore
