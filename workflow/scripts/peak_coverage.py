#!/usr/bin/env python
# Created on: Apr 9, 2024 at 12:04:43 PM
__author__ = "Michael Cuoco"

import logging, warnings

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import pandas as pd
import pyranges as pr

CHROMOSOMES = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]

cell = pd.read_parquet(snakemake.input.cell).query("max_mapq >= 30")  # type: ignore
cell["Chromosome"] = cell["Chromosome"].astype("category")
cell["Chromosome"] = cell["Chromosome"].cat.set_categories(CHROMOSOMES)
bulk = pd.read_parquet(snakemake.input.bulk)  # type: ignore


def count_overlaps(chr, group):
    """
    Count the number of overlaps between the group and the bulk data
    """

    warnings.filterwarnings("ignore", category=FutureWarning)
    my_bulk = bulk.query("Chromosome == @chr")
    brdf = pr.PyRanges(my_bulk[my_bulk["l1hs"]])
    bnrdf = pr.PyRanges(my_bulk[my_bulk["megane_gaussian"]])
    df = pr.PyRanges(group)

    return pd.DataFrame(
        {
            "l1hs_covered": df.count_overlaps(brdf).df["NumberOverlaps"].sum()
            if len(df) > 0
            else 0,
            "l1hs_missing": len(brdf)
            - df.count_overlaps(brdf).df["NumberOverlaps"].sum()
            if len(df) > 0
            else len(brdf),
            "knrgl_covered": df.count_overlaps(bnrdf).df["NumberOverlaps"].sum()
            if len(df) > 0
            else 0,
            "knrgl_missing": len(bnrdf)
            - df.count_overlaps(bnrdf).df["NumberOverlaps"].sum()
            if len(df) > 0
            else len(bnrdf),
        },
        index=[0],
    )


cell = pd.concat(
    [
        count_overlaps(chr, group)
        for chr, group in cell.groupby(["Chromosome"], observed=False)
    ]
)
cell.to_parquet(snakemake.output[0])
