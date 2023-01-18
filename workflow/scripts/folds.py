#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

# TODO: make function and object names more concise

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pickle
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import LabelEncoder
from imblearn.under_sampling import RandomUnderSampler


sys.stderr = open(snakemake.log[0], "w")

donors = [pd.read_pickle(fn).reset_index() for fn in snakemake.input.samples]

# concatenate all cells into a single table
df = (
    pd.concat(donors)
    .sort_values(["chrom", "start", "end", "cell_id", "donor_id"])
    .set_index(["chrom", "start", "end", "cell_id", "donor_id"])
)
df = df[df["all_reads.count"] >= snakemake.params.min_reads]

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
df = df.fillna(0).reset_index()
df[features] = np.minimum(
    df[features], 4e9
)  # take minimum of features and 4e9 to avoid overflow error

# make labels, using LabelEncoder() to convert strings to integers
le = LabelEncoder()
le = le.fit(list(set(df["label"])))

# save label encoder
with open(snakemake.output.label_encoder, "wb") as f:
    pickle.dump(le, f)

# get donor_id for group split
groups = df[snakemake.params.split_by]

# use StratifiedGroupKFold to preserve class balance and group balance
sgkf = StratifiedGroupKFold(n_splits=snakemake.params.num_folds)

# store labels of each fold to plot class balance
features_dict, labels_dict = {}, {}
folds = []

for fold, (train_index, test_index) in enumerate(
    sgkf.split(df[features], df["label"], groups=groups)
):
    train, test = df.iloc[train_index], df.iloc[test_index]

    if snakemake.params.downsample_train:
        train, _ = RandomUnderSampler(random_state=42).fit_resample(
            train, train["label"]
        )

    if snakemake.params.downsample_test:
        test, _ = RandomUnderSampler(random_state=42).fit_resample(test, test["label"])

    # ensure all classes are represented in the splits
    for d in [train["label"], test["label"]]:
        assert len(set(d)) == len(
            set(df["label"])
        ), "not all classes represented in split"

    features_dict[fold] = {"train": train[features], "test": test[features]}
    labels_dict[fold] = {
        "train": le.transform(train["label"]),
        "test": le.transform(test["label"]),
    }

    # append to metrics
    train["fold"] = fold + 1
    train["stage"] = "train"
    folds.append(train)
    test["fold"] = fold + 1
    test["stage"] = "test"
    folds.append(test)

# save
with open(snakemake.output.features, "wb") as f:
    pickle.dump(features_dict, f)

with open(snakemake.output.labels, "wb") as f:
    pickle.dump(labels_dict, f)

pd.concat(folds).to_pickle(snakemake.output.folds)

sys.stderr.close()
