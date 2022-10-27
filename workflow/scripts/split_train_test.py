#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

# TODO: make function and object names more concise

import sys
import pandas as pd
import numpy as np
import re
import functools
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn import metrics
import pdb


@functools.lru_cache()
def read_reference_l1():
    df = pd.read_csv(snakemake.input.ref_l1[0], index_col=[0, 1, 2])
    return df


@functools.lru_cache()
def read_non_ref_db():
    df = pd.read_csv(snakemake.input.non_ref_l1[0], index_col=[0, 1, 2])
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

    pred_df = pd.DataFrame()  # create df to predictions for each fold

    for fold, (train_index, test_index) in enumerate(sgkf.split(X, y, groups=chroms)):
        print(f"Fold {fold+1} of {num_folds}", file=sys.stderr)
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y[train_index], y[test_index]
        print(
            f"""
            training on chromosomes: {set(X_train.index.get_level_values('chrom'))}\n
            Train_x Shape: {X_train.shape}\n
            Train_y Shape: {y_train.shape}\n
            testing on chromosomes {set(X_test.index.get_level_values('chrom'))}\n
            Test_x Shape: {X_test.shape}\n
            Test_y Shape: {y_test.shape}\n
            """,
            file=sys.stderr,
        )

        # TODO: make classifier customizable by passing it as an argument
        cla = RandomForestClassifier(
            bootstrap=True,
            n_estimators=100,  # TODO: test this parameter
            oob_score=True,
            n_jobs=-1,
        )
        cla = cla.fit(X_train, y_train)  # train

        # train predictions
        train_df = X_train.copy()
        train_df["class"] = le.inverse_transform(y_train)
        train_df["pred"] = le.inverse_transform(
            cla.predict(
                X_train,
            )
        )
        train_df["phase"] = "train"
        pred_df = pred_df.append(train_df)

        # test predictions
        test_df = X_test.copy()
        test_df["class"] = le.inverse_transform(y_test)
        test_df["pred"] = le.inverse_transform(cla.predict(X_test))
        test_df["phase"] = "test"
        pred_df = pred_df.append(test_df)

        # generate metrics
        # TODO: generate metrics for training and testing
        # TODO: report feature importance
        # print(metrics.confusion_matrix(y_test, y_pred))
        # report = metrics.classification_report(
        #     y_test, y_pred, target_names=le.classes_, output_dict=True
        # )
        # del report["accuracy"]
        # for key in report.keys():
        #     report[key]["fold"] = fold + 1
        # metrics_df = metrics_df.append(
        #     pd.DataFrame(report).T.reset_index().rename(columns={"index": "class"})
        # )

    # save metrics to file
    pred_df.to_csv(snakemake.output.pred, index=True, sep="\t")


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    main(
        snakemake.input.samples, snakemake.params.num_folds, snakemake.params.min_reads
    )
    sys.stderr.close()
