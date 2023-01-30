#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys, pickle, os
from gzip import GzipFile
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold
from imblearn.under_sampling import RandomUnderSampler


def read_prep_data(files):

    # read in data
    df = pd.concat([pd.read_parquet(fn) for fn in files])

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
                "label",
                "build",
                "db",
            ]
        )
    ]

    # replace NA values with 0
    df[features] = df[features].fillna(0)

    # take minimum of features and 4e9 to avoid overflow error
    df[features] = np.minimum(df[features], 4e9)

    return df, features


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    # read and prepare data
    df, features = read_prep_data(snakemake.input)
    df = (
        df[df["all_reads.count"] >= snakemake.params.min_reads]
        .reset_index()
        .drop("index", axis=1)
    )

    # setup folds strategy
    kf = StratifiedGroupKFold(
        n_splits=snakemake.params.num_folds, shuffle=True, random_state=42
    )

    # create folds
    groups = df[snakemake.params.split_by]

    folds = {}
    for f, (train_index, test_index) in enumerate(
        kf.split(df, df["label"], groups=groups)
    ):
        folds[f] = {}

        # make train set
        folds[f]["train"] = df.loc[train_index, :]

        # make test set
        if snakemake.params.test_sampling_strategy:
            folds[f]["test"], _ = RandomUnderSampler(
                sampling_strategy=snakemake.params.test_sampling_strategy,
                random_state=42,
            ).fit_resample(df.loc[test_index, :], df.loc[test_index, "label"])
        else:
            folds[f]["test"] = df.loc[test_index, :]

    # save
    with open(snakemake.output.features, "w") as f:
        f.write("\n".join(features))

    with GzipFile(snakemake.output.folds, "wb") as f:
        pickle.dump(folds, f, protocol=-1)

    sys.stderr.close()
