#!/usr/bin/env python
# Created on: Apr 17, 2024 at 5:11:29 PM
__author__ = "Michael Cuoco"

import logging

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import pandas as pd
import pyranges as pr


def read_rmsk_bed(file: str):
    "read rmsk bed file into dataframe"

    coord_conv = lambda x: int(x.rstrip(")").lstrip("("))

    return pd.read_csv(
        file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        names=[
            "Chromosome",
            "Start",
            "End",
            "repName",
            "Score",
            "Strand",
            "milliDiv",
            "milliDel",
            "milliIns",
            "genoLeft",
            "repClassFamily",
            "repStart",
            "repEnd",
            "repLeft",
        ],
        converters={
            "repStart": coord_conv,
            "repLeft": coord_conv,
            "genoLeft": coord_conv,
        },
    )


def merge(df):
    return pd.Series(
        {
            "Chromosome": df["Chromosome"].iloc[0],
            "Start": df["Start"].min(),
            "End": df["End"].max(),
            "repStart": df["repStart"][df["repStart"].idxmin()],
            "repEnd": df["repEnd"][df["repStart"].idxmax()],
        }
    )


rmsk = read_rmsk_bed(snakemake.input[0])
rmsk = (
    pr.PyRanges(rmsk)
    .extend(200)
    .cluster()
    .df.groupby(["Cluster", "repName"], observed=True)
    .apply(merge)
    .query("repEnd > 860")
    .reset_index("repName")
    .rename(columns={"repName": "Name"})
    .reset_index(drop=True)
)

for name in ["l1hs", "l1pa2", "l1pa3", "l1pa4", "l1pa5", "l1pa6"]:
    g = name.upper() + "_3end"
    df = rmsk.query("Name == @g")
    print(f"writing {name}")
    pr.PyRanges(df).to_bed(snakemake.output[name])  # type: ignore
