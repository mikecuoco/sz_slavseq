#!/usr/bin/env python
# Created on: Aug 7, 2023 at 10:37:48 AM
__author__ = "Michael Cuoco"

import numpy as np
import pandas as pd


class Evaluator:
    def __init__(
        self, data: pd.DataFrame, proba: np.ndarray[float], label_col: str
    ) -> None:
        """
        Parameters
        ----------
        data : pd.DataFrame, data to evaluate model performance on
        proba : np.ndarray, predicted probabilities
        label_col : str, column with labels
        """
        # check inputs
        for c in [label_col, f"{label_col}_id", "cell_id"]:
            assert c in data.columns, f"test must have column {c}"
        assert (
            data.shape[0] == proba.shape[0]
        ), "test and test_proba must have the same number of rows"
        assert data[label_col].dtype == bool, f"{label_col} must be boolean"

        self.data = data[[label_col, f"{label_col}_id", "cell_id"]].copy()
        self.proba = proba
        self.label_col = label_col

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

        # get predictions for each threshold, aggregate to locus level
        locus_pred = (
            self.data[self.data[f"{self.label_col}_id"] > 0]
            .groupby([f"{self.label_col}_id", "cell_id"])
            .agg({c: "max" for c in pred_cols})
            .reset_index()
        )

        # compute recall
        tp = locus_pred[pred_cols].sum(axis=0)  # sum across loci for each threshold
        fn = locus_pred.shape[0] - tp
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

    def all_metrics(
        self, extra_pos_windows: int = 0, extra_pos_loci: int = 0
    ) -> dict[str, np.ndarray]:
        """
        Compute all metrics for each threshold

        Parameters
        ----------
        extra_pos_windows : int, optional, extra false negative windows to add to the denominator
        extra_pos_loci : int, optional, extra false negative loci to add to the denominator
        """

        self.clean_pred_cols()
        precision, recall, adjusted_recall = self.precision_recall_curve(
            self.pred_cols, extra_fn=extra_pos_windows
        )
        # null_precision, null_recall, null_adjusted_recall = self.precision_recall_curve(self.null_pred_cols, extra_fn=extra_pos_windows)
        locus_recall, adjusted_locus_recall = self.locus_recall(
            self.pred_cols, extra_fn=extra_pos_loci
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
