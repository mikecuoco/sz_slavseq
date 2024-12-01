#!/usr/bin/env python
# Created on: Apr 17, 2024 at 5:11:29 PM
__author__ = "Michael Cuoco"

import logging


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
            "Strand": df["Strand"].unique()[0]
            if len(df["Strand"].unique()) == 1
            else "+/-",
            "repStart": df["repStart"][df["repStart"].idxmin()],
            "repEnd": df["repEnd"][df["repStart"].idxmax()],
        }
    )


def adjust_line1(row):
    "adjust to the last 200 bases of line"
    if row["Strand"] == "+":
        row["Start"] = row["End"] - 300
        row["End"] = row["End"] + 300
    elif row["Strand"] == "-":
        row["End"] = row["Start"] + 300
        row["Start"] = row["Start"] - 300
    else:
        logger.warning(f"Strand is not unique for {row}")
    return row


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )

    rmsk = read_rmsk_bed(snakemake.input[0])
    rmsk = (
        pr.PyRanges(rmsk)
        .extend(200)  # extend for clustering purposes
        .cluster()
        .df.groupby(["Cluster", "repName"], observed=True)
        .apply(merge)
        .reset_index("repName")
        .rename(columns={"repName": "Name"})
        .reset_index(drop=True)
    )

    reps = {
        "l1hs": "L1HS_3end",
        "l1pa2": "L1PA2_3end",
        "l1pa3": "L1PA3_3end",
        "l1pa4": "L1PA4_3end",
        "l1pa5": "L1PA5_3end",
        "l1pa6": "L1PA6_3end",
        "polyA": "(A)n",
        "polyT": "(T)n",
    }

    for i, name in reps.items():
        df = rmsk.query("Name == @name")
        if "_3end" in name:
            df = df.query("repEnd > 860").apply(adjust_line1, axis=1)

        logger.info(f"writing {i} to {snakemake.output[i]}")  # type: ignore
        pr.PyRanges(df).to_bed(snakemake.output[i])  # type: ignore
