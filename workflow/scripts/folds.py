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
import pdb


def classes_per_cell(df, outfile):
    plot_df = (
        df.value_counts(["label", "donor_id", "cell_id"])
        .to_frame("count")
        .reset_index()
    )
    fig = sns.boxplot(plot_df, x="count", y="donor_id", hue="label")
    fig.set_xscale("log")
    sns.move_legend(fig, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig(outfile, format="png", dpi=300, bbox_inches="tight")
    plt.clf()


def classes_per_donor(df, outfile):
    plot_df = df.value_counts(["label", "donor_id"]).to_frame("count").reset_index()
    fig = sns.barplot(plot_df, x="count", y="donor_id", hue="label")
    fig.set_xscale("log")
    sns.move_legend(fig, "upper left", bbox_to_anchor=(1, 1))
    plt.savefig(outfile, format="png", dpi=300, bbox_inches="tight")
    plt.clf()


def features_per_class(df, outfile):
    plot_df = df.reset_index().melt(
        id_vars=[
            "chrom",
            "start",
            "end",
            "donor_id",
            "cell_id",
            "label",
            "build",
            "db",
        ],
        var_name="feature",
    )
    fig = sns.FacetGrid(plot_df, col="feature", col_wrap=3, sharex=False)
    fig.map_dataframe(sns.boxplot, x="value", y="label", hue="label", fliersize=0)
    plt.savefig(outfile, format="png", dpi=300)
    plt.clf()


def logplot(**kwargs):
    data = kwargs.pop("data")
    ax = sns.barplot(data, **kwargs)
    ax.set_xscale("log")


def classes_per_donor_per_fold(df, outfile):
    plot_df = (
        df.value_counts(["label", "donor_id", "fold", "stage"])
        .to_frame("count")
        .reset_index()
    )
    fig = sns.FacetGrid(plot_df, col="fold", row="stage", sharey=False)
    fig.map_dataframe(logplot, x="count", y="donor_id", hue="label")
    fig.add_legend()
    plt.savefig(outfile, format="png", dpi=300)
    plt.clf()


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    donors = [pd.read_pickle(fn).reset_index() for fn in snakemake.input.samples]

    # concatenate all cells into a single table
    df = (
        pd.concat(donors)
        .sort_values(["chrom", "start", "end", "cell_id", "donor_id"])
        .set_index(["chrom", "start", "end", "cell_id", "donor_id"])
    )
    df = df[df["all_reads.count"] >= snakemake.params.min_reads]

    # metrics
    classes_per_cell(df, snakemake.output.classes_per_cell)
    classes_per_donor(df, snakemake.output.classes_per_donor)
    features_per_class(df, snakemake.output.features_per_class)

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
    X = df.fillna(0).reset_index()
    X[features] = np.minimum(
        X[features], 4e9
    )  # take minimum of features and 4e9 to avoid overflow error

    # make labels, using LabelEncoder() to convert strings to integers
    Y = df["label"]
    le = LabelEncoder()
    le = le.fit(list(set(Y)))
    y = le.transform(Y)

    # save label encoder
    with open(snakemake.output.label_encoder, "wb") as f:
        pickle.dump(le, f)

    # get donor_id for group split
    groups = df.index.get_level_values("donor_id")

    # use StratifiedGroupKFold to preserve class balance and group balance
    sgkf = StratifiedGroupKFold(n_splits=snakemake.params.num_folds)

    # store labels of each fold to plot class balance
    features_dict = {}
    labels_dict = {}
    folds_metrics = []

    for fold, (train_index, test_index) in enumerate(
        sgkf.split(X, y, groups=groups.to_list())
    ):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y[train_index], y[test_index]

        X_train, y_train = RandomUnderSampler(random_state=42).fit_resample(
            X_train, y_train
        )

        # ensure all classes are represented in the splits
        for d in [y_train, y_test]:
            assert len(set(d)) == len(set(y)), "not all classes represented in split"

        features_dict[fold] = {"train": X_train[features], "test": X_test[features]}
        labels_dict[fold] = {"train": y_train, "test": y_test}

        # append to metrics
        X_train[["fold"]] = fold + 1
        X_train[["stage"]] = "train"
        folds_metrics.append(X_train)
        X_test[["fold"]] = fold + 1
        X_test[["stage"]] = "test"
        folds_metrics.append(X_test)

    # save folds
    with open(snakemake.output.features, "wb") as f:
        pickle.dump(features_dict, f)

    with open(snakemake.output.labels, "wb") as f:
        pickle.dump(labels_dict, f)

    classes_per_donor_per_fold(
        pd.concat(folds_metrics), snakemake.output.classes_per_donor_per_fold
    )

    sys.stderr.close()
