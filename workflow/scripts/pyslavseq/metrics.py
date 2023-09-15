#!/usr/bin/env python
# Created on: Aug 7, 2023 at 10:37:48 AM
__author__ = "Michael Cuoco"

# configure logging
import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import pyranges as pr
from .preprocessing import label
from .utilities import count_loci
import seaborn as sns
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

# TODO: add null predictions
# TODO: how to import .model_selection Model without circular import?
class Evaluator:
    def __init__(
        self,
        data,
        features,
        label_col,
        estimator,
        donor_knrgls: dict[str, pr.PyRanges],
    ) -> None:
        """
        Parameters
        ----------
        data : pd.DataFrame, data to evaluate
        features : list[str], list of feature columns
        label_col : str, name of label column
        estimator : estimator to evaluate
        donor_knrgls: dict[str, pr.PyRanges], dictionary of donor knrgls
        """

        # set attributes
        self.label_col = label_col
        self.donor_knrgls = donor_knrgls
        logger.info(f"Making predictions on {data.shape[0]} windows")
        self.proba = estimator.predict_proba(data[features])[:, 1]  # get probabilities

        # subset data to save memory
        logger.info("Loading data")
        self.data = data.loc[
            :,
            [
                "Chromosome",
                "Start",
                "End",
                label_col,
                "cell_id",
                "tissue_id",
                "donor_id",
            ],
        ].copy()

        # convert to categorical for faster groupby
        for c in ["donor_id", "cell_id"]:
            self.data[c] = self.data[c].astype("category")
            self.data[c] = self.data[c].cat.remove_unused_categories()
        self.data.reset_index(drop=True, inplace=True)

        # make thresholds
        self.thresholds = np.linspace(0, 1, 20)

        # make predictions
        self.pred_cols = np.array([f"pred_{t}" for t in range(len(self.thresholds))])
        self.data.loc[:, self.pred_cols] = self.proba[:, np.newaxis] > self.thresholds

        # make random predictions
        # self.null_pred_cols = np.array(["null_pred_{}".format(i) for i in range(len(self.thresholds))])
        # for t in self.thresholds:
        # 	pred_ind = np.random.choice(self.data[self.label_col].values, size=int(self.data.shape[0] * t))
        # 	self.data[f"null_pred_{t}"] = np.isin(self.data[self.label_col].values, pred_ind)

        tp = (
            self.data[self.pred_cols].to_numpy()
            & self.data[self.label_col].to_numpy()[:, np.newaxis]
        )
        fp = (
            self.data[self.pred_cols].to_numpy()
            & ~self.data[self.label_col].to_numpy()[:, np.newaxis]
        )
        fn = (
            ~self.data[self.pred_cols].to_numpy()
            & self.data[self.label_col].to_numpy()[:, np.newaxis]
        )

        # check where both fp and tp are zero and remove
        zero = (fp.sum(axis=0) == 0) & (tp.sum(axis=0) == 0)
        logger.info(
            f"Removing {zero.sum()} thresholds with no true positives or false positives"
        )
        tp = tp[:, ~zero]
        fp = fp[:, ~zero]
        fn = fn[:, ~zero]
        self.pred_cols = self.pred_cols[~zero]
        self.thresholds = self.thresholds[~zero]

        return None

    def locus_precision_per_cell(self, n_jobs: int) -> pd.DataFrame:
        """
        Cluster windows above a certiain threshold in each cell
        """

        def compute_precision(p, t, data, label_col):
            res = {
                "tp": [],
                "fp": [],
                "precision": [],
                "cell_id": [],
                "threshold": [],
            }
            for c, tdf in data[data[p]].groupby("cell_id"):
                if tdf.empty:
                    continue
                cdf = pr.PyRanges(tdf).cluster().df  # cluster overlapping windows
                cdf.Cluster = cdf.Cluster.astype(
                    "category"
                )  # convert Cluster to categorical
                tp = cdf.groupby("Cluster")[label_col].any().sum()

                # append results
                res["tp"].append(tp)
                res["fp"].append(cdf.Cluster.nunique() - tp)
                res["precision"].append(tp / cdf.Cluster.nunique())
                res["cell_id"].append(c)
                res["threshold"].append(t)

            return res

        res = Parallel(n_jobs=n_jobs, verbose=2)(
            delayed(compute_precision)(p, t, self.data, self.label_col)
            for p, t in zip(self.pred_cols, self.thresholds)
        )
        res = pd.DataFrame(res).explode(
            ["tp", "fp", "precision", "cell_id", "threshold"]
        )
        res["fp"] = res["fp"].astype(int)
        res["tp"] = res["tp"].astype(int)
        res["precision"] = res["precision"].astype(float)

        return res

    def locus_recall_per_cell(self, n_jobs: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute the recall and adjusted recall at the locus level for each threshold
        """

        def compute_recall(d, df):
            df.cell_id = df.cell_id.cat.remove_unused_categories()
            df = label(
                df, self.donor_knrgls[d].df, name="label", add_id=True
            )  # give each locus a unique id
            df = (
                df.query("label_id > 0")
                .groupby(["label_id", "cell_id"])
                .agg(
                    {c: "max" for c in self.pred_cols}
                )  # is there at least one overlapping window in each cell?
                .fillna(False)
                .reset_index()
                .drop(columns=["label_id"])
                .groupby("cell_id")
                .sum()
                .reset_index()
                .melt(id_vars=["cell_id"], var_name="pred_col", value_name="tp")
            )
            # map threshold onto pred_col
            df["threshold"] = df["pred_col"].map(
                {p: t for p, t in zip(self.pred_cols, self.thresholds)}
            )
            df.drop(columns=["pred_col"], inplace=True)
            df["fn"] = len(self.donor_knrgls[d].df) - df["tp"]
            df["recall"] = df["tp"] / len(self.donor_knrgls[d].df)

            return df

        res = Parallel(n_jobs=n_jobs, verbose=2)(
            delayed(compute_recall)(d, df) for d, df in self.data.groupby("donor_id")
        )

        return pd.concat(res)

    def locus_precision_recall_per_cell(self, n_jobs: int) -> pd.DataFrame:

        precision = self.locus_precision_per_cell(n_jobs)
        recall = self.locus_recall_per_cell(n_jobs)

        return pd.merge(
            precision,
            recall,
            on=["cell_id", "threshold"],
            suffixes=("_precision", "_recall"),
        )

    def precision_recall_curve(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute the precision and recall at the window level for each threshold
        """

        tp = (
            self.data[self.pred_cols].to_numpy()
            & self.data[self.label_col].to_numpy()[:, np.newaxis]
        )
        fp = (
            self.data[self.pred_cols].to_numpy()
            & ~self.data[self.label_col].to_numpy()[:, np.newaxis]
        )
        fn = (
            ~self.data[self.pred_cols].to_numpy()
            & self.data[self.label_col].to_numpy()[:, np.newaxis]
        )

        precision = tp.sum(axis=0) / (tp.sum(axis=0) + fp.sum(axis=0))
        recall = tp.sum(axis=0) / (tp.sum(axis=0) + fn.sum(axis=0))

        return precision, recall, self.thresholds

    def all_metrics(self) -> dict[str, np.ndarray]:
        """
        Compute all metrics for each threshold

        Parameters
        ----------
        extra_pos_windows : int, optional, extra false negative windows to add to the denominator
        extra_pos_loci : int, optional, extra false negative loci to add to the denominator
        """

        logger.info(
            "Calculating precision, recall, and adjusted recall at the window level"
        )
        precision, recall, adjusted_recall = self.precision_recall_curve(
            self.pred_cols, extra_fn=self.extra_pos_windows
        )
        logger.info("Calculating recall and adjusted recall at the locus level")
        # null_precision, null_recall, null_adjusted_recall = self.precision_recall_curve(self.null_pred_cols, extra_fn=extra_pos_windows)
        locus_recall, adjusted_locus_recall = self.locus_recall(
            self.pred_cols, extra_fn=self.extra_pos_loci
        )
        # null_locus_recall, null_adjusted_locus_recall = self.locus_recall(self.null_pred_cols, extra_fn=extra_pos_loci)

        # assert all are the same length
        assert (
            len(precision)
            == len(recall)
            == len(adjusted_recall)
            == len(locus_recall)
            == len(adjusted_locus_recall)
            == len(self.pred_cols)
        ), "precision, recall, adjusted_recall, locus_recall, and adjusted_locus_recall must have the same length"

        return {
            "Chromosomes": self.data["Chromosome"].unique().tolist(),
            "Donors": self.data["donor_id"].unique().tolist(),
            "Tissues": self.data["tissue_id"].unique().tolist(),
            "Cells": self.data["cell_id"].unique().tolist(),
            "threshold": self.thresholds,
            "precision": precision,
            "recall": recall,
            "adjusted_recall": adjusted_recall,
            # "null_precision": null_precision,
            # "null_recall": null_recall,
            # "null_adjusted_recall": null_adjusted_recall,
            "locus_recall": locus_recall,
            "adjusted_locus_recall": adjusted_locus_recall,
            # 	"null_locus_recall": null_locus_recall,
            # 	"null_adjusted_locus_recall": null_adjusted_locus_recall,
        }


class Visualizer:
    def __init__(self, mdl) -> None:
        self.mdl = mdl
        self.results = mdl.out
        self.loss = mdl.clf._settings["metric"]
        return None

    def learning_curve(self) -> pd.DataFrame:

        res = []
        for k in self.results:
            if isinstance(k, int):
                res.append(
                    {
                        "fold": k,
                        "time_budget": self.results[k]["time_history"],
                        "best_valid_loss_history": self.results[k][
                            "best_valid_loss_history"
                        ],
                    }
                )
        res = pd.DataFrame(res).explode(["time_budget", "best_valid_loss_history"])

        g = sns.relplot(
            res,
            x="time_budget",
            y="best_valid_loss_history",
            hue="fold",
            kind="line",
        )
        g.set(xlabel="Time budget", ylabel=f"Best Validation Set Loss ({self.loss})")

        return res

    def tuning_curve(self, y_axis: str) -> pd.DataFrame:
        assert y_axis in [
            "time_budget",
            "valid_loss_history",
        ], f"invalid y_axis {y_axis}, must be one of ['time_budget', 'valid_loss_history']"

        res = []
        for k in self.results:
            if isinstance(k, int):
                res.append(
                    {
                        "fold": k,
                        "time_budget": self.results[k]["time_history"],
                        "valid_loss_history": self.results[k]["valid_loss_history"],
                        "config_history": self.results[k]["config_history"],
                    }
                )
        res = pd.DataFrame(res).explode(
            ["time_budget", "valid_loss_history", "config_history"]
        )
        res = pd.concat(
            [
                res.drop("config_history", axis=1),
                res["config_history"].apply(pd.Series),
            ],
            axis=1,
        )
        res = pd.concat(
            [
                res.drop("Current Hyper-parameters", axis=1),
                res["Current Hyper-parameters"].apply(pd.Series),
            ],
            axis=1,
        )
        res = res.drop(
            columns=[
                "Best Hyper-parameters",
                "Best Learner",
                "Current Learner",
                "Current Sample",
            ]
        ).melt(id_vars=["time_budget", "valid_loss_history", "fold"])

        g = sns.relplot(
            data=res,
            x="value",
            y=y_axis,
            hue="fold",
            col="variable",
            col_wrap=5,
            kind="line",
            facet_kws={"sharex": False},
        )
        g.set(xlabel="Value")
        if y_axis == "time_budget":
            g.set(ylabel="Time budget")
        elif y_axis == "valid_loss_history":
            g.set(ylabel=f"Validation Set Loss ({self.loss})")

        return res

    def feature_importances(self) -> pd.DataFrame:

        res = []
        for k in self.results:
            if isinstance(k, int):
                res.append(
                    {
                        "fold": k,
                        "features": self.results[k]["features"],
                        "feature_importances": self.results[k]["feature_importances"],
                    }
                )

        res = (
            pd.DataFrame(res)
            .explode(["features", "feature_importances"])
            .pivot_table(index="features", columns="fold", values="feature_importances")
            .astype(float)
        )

        # heatmap
        g = sns.heatmap(
            res,
            annot=True,
            cmap="Blues",
        )
        g.set(xlabel="Fold", ylabel="Feature")

        return res

    def precision_recall_curve(
        self,
        recall_col: str = "adjusted_locus_recall",
        summarize_folds: bool = True,
        ax=None,
    ) -> pd.DataFrame:

        res = []
        for k in self.results:
            if isinstance(k, int):
                for stage in ["train", "test"]:
                    for r in self.results[k][stage]:
                        if set(r["Tissues"]) == set(self.mdl.data.tissue_id.unique()):
                            res.append(
                                {
                                    "fold": k,
                                    "stage": stage,
                                    "precision": r["precision"],
                                    recall_col: r[recall_col],
                                    "threshold": r["threshold"],
                                }
                            )

        res = pd.DataFrame(res).explode(["precision", recall_col, "threshold"])
        res.rename(columns={recall_col: "recall"}, inplace=True)

        if summarize_folds:
            plot_precision_recall_curve_folds(res, ax=ax)
        else:
            g = sns.relplot(
                res,
                x="recall",
                y="precision",
                hue="fold",
                col="stage",
                kind="line",
                ax=ax,
            )
            g.set(xlim=(0, 1), ylim=(0, 1))

        return res

    def precision_recall_tissues(
        self, recall_col: str = "adjusted_locus_recall"
    ) -> pd.DataFrame:
        # TODO: use plotly to make interactive plot
        # TODO: add error bars to each point

        res = []
        for k in self.results:
            if isinstance(k, int):
                for r in self.results[k]["test"]:
                    if len(r["Tissues"]) == 1:
                        res.append(
                            {
                                "fold": k,
                                "tissue": r["Tissues"][0],
                                "precision": r["precision"],
                                recall_col: r[recall_col],
                            }
                        )

        res = (
            pd.DataFrame(res)
            .explode(["precision", recall_col])
            .groupby("tissue")
            .mean()
            .reset_index()
            .drop(columns=["fold"])
        )

        g = sns.relplot(
            res,
            x=recall_col,
            y="precision",
            kind="scatter",
        )
        g.set(xlim=(0, 1), ylim=(0, 1))

        return res


def plot_precision_recall_curve_folds(df, ax=None):

    # check columns
    for c in ["precision", "recall", "threshold", "fold", "stage"]:
        assert c in df.columns, f"{c} not in df.columns"
        if c in ["threshold", "stage"]:
            df[c] = df[c].astype("category")
            df[c] = df[c].cat.remove_unused_categories()

    # get mean and std across folds for each threshold
    logger.info("Calculating mean and std across folds for each threshold")
    df = (
        df.groupby(["threshold", "stage"])[["precision", "recall"]]
        .agg(["mean", "std"])
        .fillna(0)
    )

    # compute error bars
    logger.info("Computing error bars")
    df.columns = df.columns.map("-".join)
    for c in df.columns.str.split("-").str[0].unique():
        df[f"{c}-mean_plus_std"] = df[f"{c}-mean"] + df[f"{c}-std"]
        df[f"{c}-mean_minus_std"] = df[f"{c}-mean"] - df[f"{c}-std"]
    df.reset_index(inplace=True)

    if ax is None:
        ax = plt.gca()

    # plot with matplotlib
    for s, d in df.groupby("stage"):
        ax.plot(d[f"recall-mean"], d["precision-mean"], label=s)
        ax.fill_between(
            d[f"recall-mean"],
            d["precision-mean_minus_std"],
            d["precision-mean_plus_std"],
            alpha=0.2,
        )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("recall")
    ax.set_ylabel("precision")
    ax.legend()

    return ax
