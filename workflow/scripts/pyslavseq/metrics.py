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
import matplotlib.pyplot as plt

# TODO: add null predictions
# TODO: how to import .model_selection Model without circular import?
class Evaluator:
    def __init__(
        self,
        mdl,
        query_str: str,
        donor_knrgls: dict[str, pr.PyRanges],
    ) -> None:
        """
        Parameters
        ----------
        model: Model, model to evaluate
        query_str: str, query string to subset the data to the desired set
        donor_knrgls: dict[str, pr.PyRanges], dictionary of donor knrgls
        """

        # get logger from mdl object
        # import pdb; pdb.set_trace()
        self.data = mdl.data.query(query_str)  # subset data
        mdl._validate_inputs(self.data)  # validate inputs
        self.proba = mdl.clf.predict_proba(self.data[mdl.features])[
            :, 1
        ]  # get probabilities
        self.data = self.data.query(query_str)[
            [
                "Chromosome",
                "Start",
                "End",
                mdl.label_col,
                "cell_id",
                "tissue_id",
                "donor_id",
            ]
        ]  # remove feature columns

        # set attributes
        self.label_col = mdl.label_col
        self.donor_knrgls = donor_knrgls

        # get extra fn
        raw_pos_windows = mdl.raw_data.query(query_str)[self.label_col].sum()
        filtered_pos_windows = mdl.data.query(query_str)[self.label_col].sum()
        self.extra_pos_windows = raw_pos_windows - filtered_pos_windows
        raw_pos_loci = count_loci(mdl.raw_data.query(query_str), self.donor_knrgls)
        filtered_pos_loci = count_loci(mdl.data.query(query_str), self.donor_knrgls)
        self.extra_pos_loci = raw_pos_loci - filtered_pos_loci

        # make data
        self.thresholds = np.linspace(0.01, 0.99, 99)

        # make predictions
        self.pred_cols = np.array(
            ["pred_{}".format(i) for i in range(len(self.thresholds))]
        )
        self.data.loc[:, self.pred_cols] = self.proba[:, np.newaxis] > self.thresholds

        # make random predictions
        # self.null_pred_cols = np.array(["null_pred_{}".format(i) for i in range(len(self.thresholds))])
        # for t in self.thresholds:
        # 	pred_ind = np.random.choice(self.data[self.label_col].values, size=int(self.data.shape[0] * t))
        # 	self.data[f"null_pred_{t}"] = np.isin(self.data[self.label_col].values, pred_ind)

        return None

    def clean_pred_cols(self) -> None:

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
        tp = tp[:, ~zero]
        fp = fp[:, ~zero]
        fn = fn[:, ~zero]
        self.pred_cols = self.pred_cols[~zero]
        self.thresholds = self.thresholds[~zero]

    def locus_recall(
        self, pred_cols: str, extra_fn: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Compute the recall and adjusted recall at the locus level for each threshold

        Parameters
        ----------
        extra_fn : int, optional, extra false negative loci to add to the denominator
        """
        self.data["donor_id"] = self.data["donor_id"].astype(str)
        self.data["cell_id"] = self.data["cell_id"].astype(str)
        logger.info(
            f"Computing locus recall on {self.data.cell_id.nunique()} cells from {self.data.donor_id.astype(str).unique()} donors"
        )

        # give each locus a unique id
        data = (
            self.data.groupby("donor_id")
            .apply(
                lambda d: label(
                    d, self.donor_knrgls[d.name].df, name="label", add_id=True
                )
            )
            .reset_index(drop=True)
        )

        # get the predicted loci
        locus_pred = (
            data[data["label_id"] > 0]
            .groupby(["label_id", "cell_id"])
            .agg({c: "max" for c in pred_cols})
            .reset_index()
        )

        # compute recall
        tp = locus_pred[pred_cols].sum(axis=0)  # sum across loci for each threshold
        fn = locus_pred.shape[0] - tp
        recall = tp / (tp + fn)
        adjusted_recall = tp / (tp + fn + extra_fn)

        recall = tp / (tp + fn)
        adjusted_recall = tp / (tp + fn + extra_fn)

        return recall, adjusted_recall

    def precision_recall_curve(
        self, pred_cols: str, extra_fn: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute precision, recall, and adjusted recall at the window level for each threshold

        Parameters
        ----------
        extra_fn : int, optional, extra false negative windows to add to the denominator
        """

        tp = (
            self.data[pred_cols].to_numpy()
            & self.data[self.label_col].to_numpy()[:, np.newaxis]
        )
        fp = (
            self.data[pred_cols].to_numpy()
            & ~self.data[self.label_col].to_numpy()[:, np.newaxis]
        )
        fn = (
            ~self.data[pred_cols].to_numpy()
            & self.data[self.label_col].to_numpy()[:, np.newaxis]
        )

        precision = tp.sum(axis=0) / (tp.sum(axis=0) + fp.sum(axis=0))
        recall = tp.sum(axis=0) / (tp.sum(axis=0) + fn.sum(axis=0))
        adjusted_recall = tp.sum(axis=0) / (tp.sum(axis=0) + fn.sum(axis=0) + extra_fn)

        return precision, recall, adjusted_recall

    def all_metrics(self) -> dict[str, np.ndarray]:
        """
        Compute all metrics for each threshold

        Parameters
        ----------
        extra_pos_windows : int, optional, extra false negative windows to add to the denominator
        extra_pos_loci : int, optional, extra false negative loci to add to the denominator
        """

        self.clean_pred_cols()
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
    def __init__(self, results: dict) -> None:
        self.results = results
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
        g.set(xlabel="Time budget", ylabel="Best Validation Set Loss")

        return res

    def tuning_curve(self) -> pd.DataFrame:

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
            y="valid_loss_history",
            hue="fold",
            col="variable",
            col_wrap=5,
            kind="line",
            facet_kws={"sharex": False},
        )
        g.set(xlabel="Value", ylabel="Validation Set Loss")

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
        self, recall_col: str = "adjusted_locus_recall", summarize_folds: bool = True
    ) -> pd.DataFrame:

        res = []
        for k in self.results:
            if isinstance(k, int):
                for stage in ["train", "test"]:
                    for r in self.results[k][stage]:
                        if set(r["Tissues"]) == set(self.data.tissue_id.unique()):
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

        if summarize_folds:
            # get mean and std across folds for each threshold
            df = (
                res.groupby(["threshold", "stage"])[["precision", recall_col]]
                .agg(["mean", "std"])
                .fillna(0)
            )
            df.columns = df.columns.map("-".join)
            for c in df.columns.str.split("-").str[0].unique():
                df[f"{c}-mean_plus_std"] = df[f"{c}-mean"] + df[f"{c}-std"]
                df[f"{c}-mean_minus_std"] = df[f"{c}-mean"] - df[f"{c}-std"]
            df.reset_index(inplace=True)

            # plot with matplotlib
            for s, d in df.groupby("stage"):
                plt.plot(d[f"{recall_col}-mean"], d["precision-mean"], label=s)
                plt.fill_between(
                    d[f"{recall_col}-mean"],
                    d["precision-mean_minus_std"],
                    d["precision-mean_plus_std"],
                    alpha=0.2,
                )

            plt.xlabel(recall_col)
            plt.ylabel("precision")
            plt.legend()
            plt.show()

        else:

            g = sns.relplot(
                res,
                x=recall_col,
                y="precision",
                hue="fold",
                col="stage",
                kind="line",
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
