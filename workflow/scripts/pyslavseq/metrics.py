#!/usr/bin/env python
# Created on: Aug 7, 2023 at 10:37:48 AM
__author__ = "Michael Cuoco"

# configure logging
import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics._ranking import stable_cumsum, column_or_1d, assert_all_finite


def my_precision_recall_curve(
    y_true: np.ndarray, y_score: np.ndarray, extra_fn: int = 0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    My custom implementation of precision recall curve, adjusting for additional false negatives. For binary classification only.
    Adapted from https://github.com/scikit-learn/scikit-learn/blob/f07e0138b/sklearn/metrics/_ranking.py#L865.
    :param y_true: array-like of shape (n_samples,)
                                    True binary labels. If labels are not binary, pos_label should be explicitly given.
    :param y_score: array-like of shape (n_samples,)
    :param extra_fn: int, optional, extra false negatives to add to recall
    """
    # check inputs
    assert y_true.shape == y_score.shape, "y_true and y_score must have the same shape"
    assert np.array_equal(np.unique(y_true), [0, 1]), "y_true must be binary"

    y_true = column_or_1d(y_true)
    y_score = column_or_1d(y_score)
    assert_all_finite(y_true)
    assert_all_finite(y_score)

    # sort scores and corresponding truth values
    desc_score_indices = np.argsort(y_score, kind="mergesort")[::-1]
    y_score = y_score[desc_score_indices]
    y_true = y_true[desc_score_indices]

    # y_score typically has many tied values. Here we extract
    # the indices associated with the distinct values. We also
    # concatenate a value for the end of the curve.
    distinct_value_indices = np.where(np.diff(y_score))[0]
    threshold_idxs = np.r_[distinct_value_indices, y_true.size - 1]

    # accumulate the true positives with decreasing threshold
    tps = stable_cumsum(y_true)[threshold_idxs]

    # Initialize the result array with zeros to make sure that precision[ps == 0]
    # does not contain uninitialized values.
    precision = np.zeros_like(tps)
    recall = np.zeros_like(tps)

    # compute precision and recall
    fps = 1 + threshold_idxs - tps
    ps = tps + fps
    np.divide(tps, ps, out=precision, where=(ps != 0))

    fns = np.sum(y_true) - tps + extra_fn
    rs = tps + fns
    np.divide(tps, rs, out=recall, where=(rs != 0))

    # reverse the outputs so recall is decreasing
    sl = slice(None, None, -1)
    return (
        np.hstack((precision[sl], 1)),
        np.hstack((recall[sl], 0)),
        y_score[threshold_idxs][sl],
    )


def my_recall_score(y_true: np.ndarray, y_pred: np.ndarray, extra_fn=0) -> float:
    """
    My custom implementation of recall score, adjusting for additional false negatives. For binary classification only.
    Adapted from https://github.com/scikit-learn/scikit-learn/blob/f07e0138bfee41cd2c0a5d0251dc3fe03e6e1084/sklearn/metrics/_classification.py#L1592
    :param y_true: array-like of shape (n_samples,)
                                    True binary labels. If labels are not binary, pos_label should be explicitly given.
    :param y_pred: array-like of shape (n_samples,)
    :param extra_fn: int, optional, extra false negatives to add to recall
    """

    # check inputs
    assert y_true.shape == y_pred.shape, "y_true and y_pred must have the same shape"
    assert np.array_equal(np.unique(y_true), [0, 1]), "y_true must be binary"
    assert np.array_equal(np.unique(y_pred), [0, 1]), "y_pred must be binary"

    # compute recall
    tps = np.sum((y_true == 1) & (y_pred == 1))
    fns = y_true.sum() - tps

    return tps / (tps + fns + extra_fn)


def missing_germline(data: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the number of germline variants missing from the single-cell data for each chromosome.
    :param data: DataFrame with single-cell data.
    :param ref: DataFrame of reference intervals
    """

    for c in ["cell_id", "Chromosome"]:
        assert c in data.columns, f"Column '{c}' not in DataFrame"
    assert data.donor_id.nunique() == 1, "Data must be from a single donor"

    # groupby function
    def count_overlaps(group):

        # check group is unique
        out = {}
        for c in ["cell_id", "Chromosome", "donor_id"]:
            assert group[c].nunique() == 1, f"Group must be unique by '{c}'"
            out[c] = group[c].unique()[0]
            if out[c] in group.name:
                out.pop(c)

        chrom = group["Chromosome"].unique[0]
        df = pr.PyRanges(group)
        rdf = pr.PyRanges(ref.query("Chromosome == @chrom"))

        out["covered"] = df.count_overlaps(rdf).df["NumberOverlaps"].sum()
        out["missing"] = len(rdf) - out["covered"]
        out["fraction"] = out["covered"] / len(rdf)

        return pd.Series(out)

    # Apply the function to each group of cells
    return data.groupby(["cell_id", "Chromosome"]).apply(count_overlaps)


# def summarize_cns_peaks(df):
