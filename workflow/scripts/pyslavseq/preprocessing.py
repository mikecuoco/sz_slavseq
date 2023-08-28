#!/usr/bin/env python
# Created on: Aug 8, 2023 at 1:53:30 PM
__author__ = "Michael Cuoco"

# configure logging
import logging

logger = logging.getLogger(__name__)

from pathlib import Path
import pyarrow.parquet as pq
import numpy as np
import pandas as pd
import pyranges as pr


def read_cell_features(
    filename: str, cell_id: str, exclude_ref: bool = True
) -> pd.DataFrame:
    """
    Read in a single cell's windows and annotate with labels

    Parameters
    ----------
    filename : str, parquet file to read
    cell_id : str, cell id to label windows with
    """
    assert Path(filename).exists(), f"File {filename} does not exist"

    data = pq.read_table(filename).to_pandas()
    data["cell_id"] = cell_id

    # remove windows with ref reads
    if exclude_ref:
        data = data.loc[data["n_ref_reads"] == 0, :]

    return data


def label(
    df: pd.DataFrame, other_df: pd.DataFrame, name: str, add_id: bool = False
) -> pd.DataFrame:
    """
    Label windows in df with whether they overlap with other_df

    Parameters
    ----------
    df : pd.DataFrame, windows to label
    other_df : pd.DataFrame, windows to overlap with
    name : str, name of column to add to df
    add_id : bool, whether to add a column with the id of the overlapping window
    """
    # assert dfs are not empty
    assert df.shape[0] > 0, "df is empty!"
    assert other_df.shape[0] > 0, "other_df is empty!"

    # check inputs
    for c in ["Chromosome", "Start", "End"]:
        assert c in df.columns, f"df must have column {c}"
        assert c in other_df.columns, f"other_df must have column {c}"

    if add_id:
        other_df.loc[:, f"{name}_id"] = other_df.index.values
        other_df = other_df[[f"{name}_id", "Chromosome", "Start", "End"]]
    else:
        other_df = other_df[["Chromosome", "Start", "End"]]

    # join with pyranges to find overlaps
    overlap = (
        pr.PyRanges(df[["Chromosome", "Start", "End"]], int64=True)
        .join(
            pr.PyRanges(other_df, int64=True),
            how="left",
        )
        .df
    )

    overlap[name] = overlap["Start_b"] != -1
    overlap = overlap.drop(columns=["Start_b", "End_b"]).drop_duplicates()

    return df.join(
        overlap.set_index(["Chromosome", "Start", "End"]),
        on=["Chromosome", "Start", "End"],
        how="left",
    )
