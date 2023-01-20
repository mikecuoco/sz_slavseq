#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

# TODO: make function and object names more concise

import sys, os
import pandas as pd
import polars as pl
import pickle
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import LabelEncoder
from imblearn.under_sampling import RandomUnderSampler

os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
sys.stderr = open(snakemake.log[0], "w")

donors = [pl.read_parquet(fn) for fn in snakemake.input.samples]

# concatenate all cells into a single table
df = pl.concat(donors).sort(["chrom", "start", "end", "cell_id", "donor_id"])

# remove windows with too few reads
df = df.filter(pl.col("all_reads.count") >= snakemake.params.min_reads)

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

# replace NaN and null values with 0
df = df.fill_null(strategy="zero").fill_nan(0)

# take minimum of features and 4e9 to avoid overflow error
df = df.with_columns(pl.col(features).apply(lambda x: min(x, 4e9)))

# make labels, using LabelEncoder() to convert strings to integers
le = LabelEncoder()
le = le.fit(list(set(df["label"])))

# save label encoder
with open(snakemake.output.label_encoder, "wb") as f:
    pickle.dump(le, f)

# get variable for group split
groups = df[snakemake.params.split_by]

# use StratifiedGroupKFold to preserve class balance and group balance
sgkf = StratifiedGroupKFold(n_splits=snakemake.params.num_folds)

# store labels of each fold to plot class balance
features_dict, labels_dict = {}, {}
folds = []

for fold, (train_index, test_index) in enumerate(
    sgkf.split(df[features], df["label"], groups=groups)
):

    train, test = df[train_index, :].to_pandas(), df[test_index, :].to_pandas()

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
