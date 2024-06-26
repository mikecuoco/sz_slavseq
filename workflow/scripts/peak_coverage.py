#!/usr/bin/env python
# Created on: Apr 9, 2024 at 12:04:43 PM
__author__ = "Michael Cuoco"

import logging, warnings

warnings.filterwarnings("ignore", category=FutureWarning)  # silence pyranges warning

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import pandas as pd
import pyranges as pr  # type: ignore
from snakemake.io import Namedlist

CHROMOSOMES = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]

cell = pd.read_parquet(snakemake.input.cell).query("max_mapq >= 60")  # type: ignore
cell["Chromosome"] = cell["Chromosome"].astype("category")
cell["Chromosome"] = cell["Chromosome"].cat.set_categories(CHROMOSOMES)
bulk = pd.read_parquet(snakemake.input.bulk)  # type: ignore

anno = {}
for n in [
    "l1hs",
    "megane_gaussian",
    "megane_percentile",
    "megane_breakpoints",
    "graffite",
    "xtea",
]:
    assert snakemake.input.get(n), f"Missing input for annotation {n}"
    if isinstance(snakemake.input.get(n), Namedlist):  # type: ignore
        anno[n] = pr.read_bed(snakemake.input[n][0]).df  # type: ignore
    else:
        anno[n] = pr.read_bed(snakemake.input[n]).df  # type: ignore


def count_overlaps(chrom, group):
    """
    Count the number of overlaps between the group and the bulk data
    """
    chrom = chrom[0]  # for some reason chrom comes in as a tuple

    # check inputs
    if group["Chromosome"].nunique() > 1:
        logger.warning(f"Group contains chromosomes other than {chrom}")
        import pdb

        pdb.set_trace()

    # get the cell and bulk data for this chromosome
    my_cell = pr.PyRanges(group)
    my_bulk = bulk.query("Chromosome == @chrom")
    my_anno = {n: pr.PyRanges(d[d["Chromosome"] == chrom]) for n, d in anno.items()}

    # get coverage for the bulk data
    res = {}
    bdf = pr.PyRanges(my_bulk)
    res["total_bulk"] = len(bdf.df)
    res["bulk_covered"] = (
        my_cell.count_overlaps(bdf).df["NumberOverlaps"].sum() if len(my_cell.df) else 0
    )

    # get coverage for each annotation by itself and of the bulk data
    for n in anno.keys():
        res[f"total_{n}"] = len(my_anno[n].df)
        res[f"covered_{n}"] = (
            my_cell.count_overlaps(my_anno[n]).df["NumberOverlaps"].sum()
            if len(my_cell.df)
            else 0
        )
        badf = pr.PyRanges(my_bulk[my_bulk[n]])
        res[f"total_bulk_{n}"] = len(badf.df)
        res[f"covered_bulk_{n}"] = (
            my_cell.count_overlaps(badf).df["NumberOverlaps"].sum()
            if len(my_cell.df)
            else 0
        )

    return pd.DataFrame(
        index=[chrom],
        data=res,
    )


cell = pd.concat(
    [
        count_overlaps(chrom, group)
        for chrom, group in cell.groupby(["Chromosome"], observed=False)
    ]
)

cell.to_parquet(snakemake.output[0])  # type: ignore
