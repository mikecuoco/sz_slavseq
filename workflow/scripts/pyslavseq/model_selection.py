#!/usr/bin/env python
# Created on: Aug 7, 2023 at 10:37:48 AM
__author__ = "Michael Cuoco"

# configure logging
import logging

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

import json, warnings, tempfile
from typing import Union
from pathlib import Path
import numpy as np
import pandas as pd
import pyranges as pr
from sklearn.model_selection import KFold, StratifiedGroupKFold
from flaml import AutoML
from flaml.automl.data import get_output_from_log
from .metrics import Evaluator
from .utilities import count_loci


class SampleChrSplitter(KFold):
    "Class to split data by sample and chromosome, stratified to preserve class balance"

    def __init__(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        sample_col: str,
        n_chr_splits: int = 2,
        n_sample_splits: int = 2,
        shuffle: bool = False,
        random_state: int | np.random.RandomState | None = None,
    ) -> None:
        for c in [sample_col, "Chromosome"]:
            assert c in X.columns, f"X must have column {c}"

        super().__init__(
            n_splits=int(n_chr_splits * n_sample_splits),
            shuffle=shuffle,
            random_state=random_state,
        )

        self.sample_col = sample_col
        self.sample_array = X[sample_col].values
        self.chr_array = X["Chromosome"].values
        self.X = X
        self.y = y
        self.n_chr_splits = n_chr_splits
        self.n_sample_splits = n_sample_splits
        self.shuffle = shuffle

        return None

    def split(self, X):
        assert (
            X.shape[0] == self.y.shape[0]
        ), "X and y must have the same number of rows"

        sgkf = StratifiedGroupKFold(n_splits=self.n_chr_splits)
        chr_split = sgkf.split(X, self.y, groups=self.chr_array)
        sgkf = StratifiedGroupKFold(n_splits=self.n_sample_splits)
        s_split = sgkf.split(X, self.y, groups=self.sample_array)

        for chr_train, chr_test in chr_split:
            for sample_train, sample_test in s_split:
                train = np.intersect1d(sample_train, chr_train)
                test = np.intersect1d(sample_test, chr_test)

                # shuffle if requested
                if self.shuffle:
                    train = train.sample(frac=1, random_state=self.random_state)
                    test = test.sample(frac=1, random_state=self.random_state)

                yield train, test

    def get_n_splits(self):
        return self.n_splits


class NpEncoder(json.JSONEncoder):
    "Encode pandas series and numpy arrays for json output"

    def default(self, obj):
        # convert categorical columns to strings
        if isinstance(obj, pd.Categorical):
            return obj.astype(str).tolist()
        if isinstance(obj, pd.Series):
            return obj.to_numpy().tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


class Model:
    def _validate_inputs(
        self,
        data: pd.DataFrame,
    ):
        "Validate inputs to Model class"

        # check inputs
        for c in [
            self.label_col,
            "cell_id",
            "donor_id",
            "tissue_id",
            "Chromosome",
            "Start",
            "End",
        ]:
            assert c in data.columns, f"data must have column {c}"
        for f in self.features:
            assert f in data.columns, f"data must have feature column {f}"
        assert data[self.label_col].dtype == bool, f"{self.label_col} must be boolean"

        # check we have the right donor IDs
        assert set(data.donor_id.unique()).issubset(
            set(self.donor_knrgls.keys())
        ), "Some donors in data do not have corresponding knrgls"

        # check all donor_knrgls are in pyranges format
        for f in self.donor_knrgls.values():
            assert isinstance(f, pr.PyRanges), f"{f} must be a pyranges object"

    def __init__(
        self,
        clf: AutoML,
        data: pd.DataFrame,
        features: list,
        label_col: str,
        donor_knrgls: dict[str : pr.PyRanges],
        query_str: str,
        outfile: Union[str, None] = None,
    ) -> None:
        """
        Parameters
        ----------
        data : pd.DataFrame, data to fit model on
        features : list, feature columns to use for model
        label_col : str, column with labels
        donor_knrgls : dict, {donor_id: knrgl bed file in pyranges format}
        query_str : str, query string to filter data
        outfile : str, file to save results to
        """

        # save attributes
        self.clf = clf
        self.query_str = query_str
        self.features = features
        self.label_col = label_col
        self.donor_knrgls = donor_knrgls
        self.raw_data = data
        self.outfile = outfile

        # validate raw data
        self._validate_inputs(data)

        # filter data
        logger.info(f"Keeping windows with {query_str}")
        self.data = data.query(query_str).copy()  # TODO: try removing this
        # initialize proba columns
        self.data["train_proba"] = np.nan
        self.data["test_proba"] = np.nan

        # validate filtered data
        self._validate_inputs(self.data)

        # check outfile
        if self.outfile:
            assert not Path(
                self.outfile
            ).exists(), f"{self.outfile} already exists! Please choose a different name or delete the file if you want to overwrite it"

        # save info to log
        self.out = {
            "Input data": {
                "Number of total windows": len(self.raw_data),
                "Filter": query_str,
                "Number of tfiltered windows": len(self.data),
                "Features": self.features,
                "Class label": self.label_col,
                "Number of cells": self.data.cell_id.nunique(),
                "Cells": self.data.cell_id.unique(),
                "Number of tissues": self.data.tissue_id.nunique(),
                "Tissues": self.data.tissue_id.unique(),
                "Number of donors": self.data.donor_id.nunique(),
                "Donors": self.data.donor_id.unique(),
            },
            "Model settings": clf._settings,
        }

        return None

    def tune(
        self, train_idx: np.ndarray, n_hyper_splits: int, log_file_name: str
    ) -> dict[str, np.ndarray]:
        """
        Tune a model's hyperparameters, refit on training data, test it on held out data and compute metrics

        Parameters
        ----------
        train_idx : np.ndarray, indices of training set
        n_hyper_splits : int, number of splits for hyperparameter tuning
        log_file_name : str, name of log file to save results to
        """

        # define train sets
        train = self.data.iloc[train_idx, :].reset_index()

        logger.info(
            f"""
		Tuning model with {len(self.features)} features on {train.shape[0]} windows
		{train[self.label_col].sum()} positive windows ({count_loci(train, self.donor_knrgls)} loci)
		{len(train) - train[self.label_col].sum()} negative windows
		{train.Chromosome.nunique()} Chromosomes: {train.Chromosome.astype(str).unique()}
		{train.donor_id.nunique()} Donors: {train.donor_id.astype(str).unique()}
		{train.tissue_id.nunique()} Tissues: {train.tissue_id.astype(str).unique()}
		{train.cell_id.nunique()} cells
		"""
        )

        # suppress UserWarning: `use_label_encoder` is deprecated in 1.7.0.
        warnings.filterwarnings("ignore", category=UserWarning)
        self.clf.fit(
            train[self.features],
            train[self.label_col],
            split_type=StratifiedGroupKFold(n_splits=n_hyper_splits),
            groups=train.Chromosome,
            log_file_name=log_file_name,
        )
        logger.info("Done tuning")

    def pred_eval(
        self, train_idx: np.ndarray, test_idx: np.ndarray
    ) -> dict[str, np.ndarray]:

        logger.info("Evaluting model")

        # make predictions for train and test sets
        result = []
        for name, idx in {"train": train_idx, "test": test_idx}.items():
            df = self.data.iloc[idx, :].reset_index()
            # get all donor knrgls that overlap this split
            chrs = df["Chromosome"].unique().tolist()
            donor_knrgls = {
                k: df[df.Chromosome.isin(chrs)] for k, df in self.donor_knrgls.items()
            }

            # get metrics
            # for each tissue
            if df.tissue_id.nunique() > 1:
                for t, d in df.groupby("tissue_id"):
                    logger.info(f"Evaluating {name} set for tissue {t}")
                    metrics = Evaluator(
                        self, f"tissue_id == '{t}' & Chromosome in {chrs}", donor_knrgls
                    ).all_metrics()
                    metrics["stage"] = name
                    metrics["n_windows"] = len(d)
                    result.append(metrics.copy())

            # for everything
            logger.info(f"Evaluating {name} set for tissues {df.tissue_id.unique()}")
            metrics = Evaluator(
                self, f"Chromosome in {chrs}", donor_knrgls
            ).all_metrics()
            metrics["stage"] = name
            metrics["n_windows"] = len(df)
            result.append(metrics.copy())

            # add predictions
            self.data.iloc[idx].loc[:, f"{name}_proba"] = self.clf.predict_proba(
                df[self.features]
            )[:, 1]

        logger.info("Done evaluating")

        return tuple(result)

    def cv(self, n_splits: int) -> None:
        logger.info("Starting cross-validation")

        # define splitter from cross-validation
        sgkf = StratifiedGroupKFold(n_splits=n_splits)
        # convert Chromosome to category for faster grouping
        self.data.loc[:, "Chromosome"] = self.data["Chromosome"].astype("category")

        for i, (train_idx, test_idx) in enumerate(
            sgkf.split(
                self.data, self.data[self.label_col], groups=self.data["Chromosome"]
            )
        ):
            logger.info(f"Fold {i+1}/{n_splits}")
            self.out[i + 1] = {}
            flog = self.out[i + 1]

            # tune model
            with tempfile.NamedTemporaryFile() as tmp:
                self.tune(train_idx, n_hyper_splits=5, log_file_name=tmp.name)

                (
                    time_history,
                    best_valid_loss_history,
                    valid_loss_history,
                    config_history,
                    metric_history,
                ) = get_output_from_log(filename=tmp.name, time_budget=1e6)

            # save results
            flog["train"], flog["test"] = [], []
            for r in self.pred_eval(train_idx, test_idx):
                if r["stage"] == "train":
                    flog["train"].append(r)
                elif r["stage"] == "test":
                    flog["test"].append(r)
                else:
                    raise ValueError(f"Cannot read result {r}")

            # save model fit info
            flog.update(
                {
                    "time_history": time_history,
                    "best_valid_loss_history": best_valid_loss_history,
                    "valid_loss_history": valid_loss_history,
                    "config_history": config_history,
                    "metric_history": metric_history,
                    "feature_importances": self.clf.feature_importances_,
                    "features": self.clf.feature_names_in_,
                    "best_estimator": self.clf.best_estimator,
                    "best_config": self.clf.best_config,
                    "total_loci_train": count_loci(
                        self.data.iloc[train_idx], self.donor_knrgls
                    ),
                }
            )
        logger.info("Done cross-validation")

        if self.outfile:
            with open(self.outfile, "w") as f:
                json.dump(self.out, f, cls=NpEncoder)
            logger.info(f"Results saved to {self.outfile}")
