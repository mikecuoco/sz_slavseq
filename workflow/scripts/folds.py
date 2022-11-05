#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

# TODO: make function and object names more concise
# TODO: make splitting deterministic

import sys
import pandas as pd
import numpy as np
import re
import functools
import pickle
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import LabelEncoder
from src.genome.windows import read_rmsk, make_l1_windows


@functools.lru_cache()
def read_reference_l1():
    df = read_rmsk(snakemake.input.ref_l1)
    df = make_l1_windows(df, snakemake.input.chromsizes, "reference_l1hs_l1pa2_6")
    return df


@functools.lru_cache()
def read_non_ref_db():
    df = pd.read_csv(snakemake.input.non_ref_l1[0], sep = "\t", header = None, names = ["chrom", "start", "end"])
    df = make_l1_windows(df, snakemake.input.chromsizes, "in_NRdb")
    return df


def read_cell_features(fn, cell_id):
    """read a cell's feature table from pickle file"""
    # TODO make the non-ref and ref-l1 column names flexible to changes
    df = (
        pd.read_pickle(fn)
        .merge(read_non_ref_db(), left_index=True, right_index=True, how="left")
        .merge(read_reference_l1(), left_index=True, right_index=True, how="left")
        .fillna({"in_NRdb": False, "reference_l1hs_l1pa2_6": False})
    )

    df["cell_id"] = cell_id

    return df


# TODO: add cell id and donor to table in feature extraction step, then remove this function
def cells_from_sample(sample_files):
    for fn in sample_files:
        m = re.search("/([^/]+?)\.pickle\.gz$", fn)
        yield fn, m.group(1)


def label(df):
    for x in df[["in_NRdb", "reference_l1hs_l1pa2_6"]].itertuples():
        if x.reference_l1hs_l1pa2_6:
            yield "RL1"
        elif x.in_NRdb:
            yield "KNRGL"
        else:
            yield "OTHER"


def main(files, num_folds, min_reads):

    # read in feature tables for each cell
    cells = []
    for fn, cell_id in cells_from_sample(files):
        cells.append(read_cell_features(fn, cell_id).reset_index())

    # concatenate all cells into a single table, remove windows below min_reads
    df = (
        pd.concat(cells)
        .sort_values(["chrom", "start", "end", "cell_id"])
        .set_index(["chrom", "start", "end", "cell_id"])
    )
    df = df[df["all_reads.count"] >= min_reads]

    # make features
    features = [
        x
        for x in df.columns
        if not any(
            x.endswith(y)
            for y in [
                "chrom",
                "start",
                "end",
                "cell_id",
                "in_NRdb",
                "reference_l1hs_l1pa2_6",
                "fold",
                "_peak_position",
                "_en_motif",
                "_te_strand",
            ]
        )
    ]
    X = df[features].fillna(0)
    X = np.minimum(X, 4e9)  # take minimum of features and  4e9 to avoid overflow

    # make labels, using LabelEncoder() to convert strings to integers
    Y = pd.Series(label(df), index=df.index)
    le = LabelEncoder()
    le = le.fit(list(set(Y)))
    y = le.fit_transform(Y)

    # get chromosome names for group split
    chroms = df.index.get_level_values("chrom").to_list()

    # use groupkfold to split by chromosome and train/test
    sgkf = StratifiedGroupKFold(n_splits=num_folds)

    for fold, (train_index, test_index) in enumerate(sgkf.split(X, y, groups=chroms)):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # save train/test splits
        X_train.to_pickle(snakemake.output.train_features[fold])
        X_test.to_pickle(snakemake.output.test_features[fold])

        with open(snakemake.output.train_labels[fold], "wb") as f:
            pickle.dump(y_train, f)

        with open(snakemake.output.test_labels[fold], "wb") as f:
            pickle.dump(y_test, f)

    # save label encoder
    with open(snakemake.output.label_encoder, "wb") as f:
        pickle.dump(le, f)

if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    main(
        snakemake.input.samples, snakemake.params.num_folds, snakemake.params.min_reads
    )
    sys.stderr.close()
