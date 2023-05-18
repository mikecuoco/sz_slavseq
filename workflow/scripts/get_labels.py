#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import functools, sys
import numpy as np
import pandas as pd
import pyranges as pr


@functools.lru_cache()
def read_rmsk(rmsk_outfile):
    """
    Read the repeatmasker output table and return locations of L1HS and L1PA2-6
    """
    # read the rmsk file
    df = pd.read_csv(
        rmsk_outfile,
        skiprows=3,
        delim_whitespace=True,
        names=["Chromosome", "Start", "End", "Strand", "repeat"],
        dtype={"Chromosome": str, "Start": int, "End": int},
        usecols=[4, 5, 6, 8, 9],
    )

    # filter for rep_names
    rep_names = [
        "L1HS_3end",
        "L1PA2_3end",
        "L1PA3_3end",
        "L1PA4_3end",
        "L1PA5_3end",
        "L1PA6_3end",
    ]

    df = df[df["repeat"].isin(rep_names)]

    # save to new dataframe
    df["Strand"] = df["Strand"].str.replace("C", "-")
    df["Start"] -= 1  # make zero-based

    # extend the repeats by 1000 bp
    df["Start"] = df.apply(
        lambda x: x["Start"] - 1000 if x["Strand"] == "-" else x["Start"], axis=1
    )
    df["End"] = df.apply(
        lambda x: x["End"] + 1000 if x["Strand"] == "+" else x["End"], axis=1
    )

    return df


@functools.lru_cache()
def read_knrgl(knrgl_bedfile):
    df = pd.read_csv(
        knrgl_bedfile,
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "End", "Strand"],
        dtype={"Chromosome": str, "Start": int, "End": int},
        usecols=[0, 1, 2, 3],
    )

    df["Start"] = df.apply(
        lambda x: x["Start"] - 1000 if x["Strand"] == "-" else x["Start"], axis=1
    )
    df["End"] = df.apply(
        lambda x: x["End"] + 1000 if x["Strand"] == "+" else x["End"], axis=1
    )

    return df


def label_windows(df: pd.DataFrame, other_df: pd.DataFrame, label: str):
    assert type(df) == pd.DataFrame, "df must be a pandas DataFrame"
    assert type(other_df) == pd.DataFrame, "other_df must be a pandas DataFrame"
    assert label not in df.columns, f"{label} already in df.columns"

    # convert to pyranges
    pr_df = pr.PyRanges(df)
    pr_other_df = pr.PyRanges(other_df)

    # get the windows that overlap with the other_df
    overlapping = pr_df.overlap(pr_other_df).df

    # set the index to the chromosome, start, and end
    # TODO: check if reads are in same orientation as repeats
    df.set_index(["Chromosome", "Start", "End"], inplace=True)
    overlapping.set_index(["Chromosome", "Start", "End"], inplace=True)

    # label the windows that overlap
    df[label] = df.index.isin(overlapping.index)

    # reset the index
    df.reset_index(inplace=True)

    return df


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    # get features all the cells from this donor
    df = pd.concat([pd.read_parquet(f) for f in snakemake.input.features])

    # read in the repeatmasker, knrgl, and blacklist files
    rmsk_df = read_rmsk(snakemake.input.rmsk)
    knrgl_df = read_knrgl(snakemake.input.knrgl)

    # label the windows for classification
    df = label_windows(df, knrgl_df, "knrgl")
    df = label_windows(df, rmsk_df, "rmsk")

    label = []
    for row in df.itertuples():
        if row.rmsk:
            label.append("RMSK")
        elif row.knrgl:
            label.append("KNRGL")
        else:
            label.append("OTHER")

    df.drop(["knrgl", "rmsk"], axis=1, inplace=True)
    df["label"] = label

    # save
    df.to_parquet(snakemake.output[0], index=False)

    sys.stderr.close()
