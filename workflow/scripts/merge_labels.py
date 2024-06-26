#!/usr/bin/env python
# Created on: Nov 17, 2023 at 6:45:57 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

from pathlib import Path
import pandas as pd


def coverage(df):
    """
    Get the coverage of the annotation in the regions
    """

    # count how many rows are not periods
    tot = len(df)
    cov = df.loc[:, chr_col].str.contains("chr").sum()
    return pd.Series({"total": tot, "covered": cov})


if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )

    logger = logging.getLogger(__name__)

    regions = pd.read_parquet(snakemake.input.regions).set_index(
        ["Chromosome", "Start", "End"]
    )

    cov = {}

    for a in map(Path, snakemake.input.annotations):  # type: ignore
        name = a.name.split(".")[1]
        # check if file is empty
        if a.stat().st_size == 0:
            regions[name] = False
            continue
        anno = pd.read_csv(a, sep="\t", header=None)
        # find columns with chromosome, start, and end
        chr_cols = anno.apply(lambda x: x.str.contains("chr"), axis=1).sum() > 0
        # get second true value, check if this doesn't exist
        if chr_cols.sum() < 2:
            logger.warning(
                f"{name} annotation does not contain any intersection with regions."
            )
            regions[name] = False
            continue

        chr_col = chr_cols[chr_cols].index[1]
        cov[name] = anno.groupby(0).apply(lambda x: coverage(x), include_groups=False)
        # get the chromosome, start, and end columns
        anno = anno.loc[:, chr_col : (chr_col + 2)]
        anno.columns = ["Chromosome", "Start", "End"]
        anno.set_index(["Chromosome", "Start", "End"], inplace=True)

        # find regions that are in the annotation
        regions[name] = regions.index.isin(anno.index)

    (
        pd.concat(cov.values(), keys=cov.keys())
        .reset_index()
        .rename(columns={"level_0": "annotation", 0: "Chromosome"})
        .to_csv(snakemake.output.coverage, sep="\t", index=False)  # type: ignore
    )

    regions.reset_index().to_parquet(snakemake.output.data)  # type: ignore
