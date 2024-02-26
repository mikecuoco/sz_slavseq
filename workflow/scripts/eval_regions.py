#!/usr/bin/env python
# Created on: Feb 14, 2024 at 9:54:43 AM
__author__ = "Michael Cuoco"

import logging

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import warnings

warnings.filterwarnings("ignore")

from pathlib import Path
import pyranges as pr
import pandas as pd


# 1. How many annotations are covered by bulk SLAVseq bam?
# 2. How many bulk-bam-covered annotations are covered by bulk SLAVseq regions?
# 3. How many sc-bam-covered-annotations are covered by single-cell SLAVseq bam?
# 4. How many bulk-bam-covered-annotations are covered by single-cell SLAVseq bam?
# 5. How many bulk SLAVseq regions are covered by single-cell SLAVseq bam?
# 6. How many bulk SLAVseq regions are covered by single-cell SLAVseq regions?


cell_regions = pd.read_parquet(snakemake.input.cell_regions)  # type: ignore

ov = {}
for f in snakemake.input.cell_coverage:  # type: ignore
    name = Path(f).name.split(".")[1]

    anno = pd.read_csv(f, sep="\t", header=None)
    anno = anno.iloc[:, list(range(3)) + list(range(-4, 0))]
    # set names
    anno.columns = [
        "Chromosome",
        "Start",
        "End",
        "Count",
        "Breadth",
        "Width",
        "Fraction",
    ]
    anno = anno.query("Count >= 5").reset_index(drop=True)

    # read regions
    regions = cell_regions[cell_regions[name]]
    ov[(name, "regions_covered_5reads")] = len(
        pr.PyRanges(cell_regions).overlap(pr.PyRanges(anno)).df
    )
    ov[(name, "labelled_regions_covered_5reads")] = len(
        pr.PyRanges(regions).overlap(pr.PyRanges(anno)).df
    )
    ov[(name, "total_covered_5reads")] = len(anno)


if "gDNA" not in snakemake.input.cell_regions:  # type: ignore
    anno = pd.read_parquet(snakemake.input.bulk_regions)  # type: ignore
    regions = cell_regions[cell_regions["bulk"]]
    ov[("bulk", "regions_bulk")] = len(
        pr.PyRanges(cell_regions).overlap(pr.PyRanges(anno)).df
    )
    ov[("bulk", "labelled_regions_bulk")] = len(
        pr.PyRanges(regions).overlap(pr.PyRanges(anno)).df
    )
    ov[("bulk", "total_bulk")] = len(anno)

df = pd.DataFrame.from_dict(ov, orient="index")
df.columns = ["Count"]
df.reset_index(inplace=True)
df[["Name", "Type"]] = pd.DataFrame(df["index"].tolist(), index=df.index)
df.drop(columns=["index"], inplace=True)
df.pivot(columns="Type", values="Count", index="Name").reset_index(inplace=True)

df.to_csv(snakemake.output[0], sep="\t", index=False)  # type: ignore
