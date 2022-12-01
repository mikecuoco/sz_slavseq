#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

# TODO: make function and object names more concise

import sys
import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import LabelEncoder
from imblearn.over_sampling import RandomOverSampler


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
    for fn in files:
        cells.append(pd.read_pickle(fn).reset_index())

    # concatenate all cells into a single table, remove windows below min_reads
    df = (
        pd.concat(cells)
        .sort_values(["chrom", "start", "end", "cell_id", "donor_id"])
        .set_index(["chrom", "start", "end", "cell_id", "donor_id"])
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
                "donor_id",
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
    X = np.minimum(X, 4e9)  # take minimum of features and 4e9 to avoid overflow error

    # make labels, using LabelEncoder() to convert strings to integers
    Y = pd.Series(label(df), index=df.index)
    le = LabelEncoder()
    le = le.fit(list(set(Y)))
    y = le.transform(Y)

    # save label encoder
    with open(snakemake.output.label_encoder, "wb") as f:
        pickle.dump(le, f)

    # get cell_id for group split
    # TODO: split by donor
    groups = df.index.get_level_values("cell_id")

    # use groupkfold to split by chromosome and train/test
    sgkf = StratifiedGroupKFold(n_splits=num_folds)

    for fold, (train_index, test_index) in enumerate(
        sgkf.split(X, y, groups=groups.to_list())
    ):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y[train_index], y[test_index]

        # ensure all classes are represented in the splits
        for d in [y_train, y_test]:
            assert len(set(d)) == len(set(y)), "not all classes represented in split"

        X_train, y_train = RandomOverSampler(random_state=42).fit_resample(
            X_train, y_train
        )
        X_test, y_test = RandomOverSampler(random_state=42).fit_resample(X_test, y_test)

        # save train/test splits
        X_train.to_pickle(snakemake.output.train_features[fold])
        X_test.to_pickle(snakemake.output.test_features[fold])

        with open(snakemake.output.train_labels[fold], "wb") as f:
            pickle.dump(y_train, f)

        with open(snakemake.output.test_labels[fold], "wb") as f:
            pickle.dump(y_test, f)


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    main(
        snakemake.input.samples, snakemake.params.num_folds, snakemake.params.min_reads
    )
    sys.stderr.close()
