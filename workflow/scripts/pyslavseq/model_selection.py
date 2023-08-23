#!/usr/bin/env python
# Created on: Aug 7, 2023 at 10:37:48 AM
__author__ = "Michael Cuoco"

import json, warnings, tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, StratifiedGroupKFold
from .metrics import Evaluator
from flaml import AutoML
from flaml.automl.data import get_output_from_log


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

        for chr_train, chr_test in StratifiedGroupKFold(
            n_splits=self.n_chr_splits
        ).split(X, self.y, groups=self.chr_array):
            for sample_train, sample_test in StratifiedGroupKFold(
                n_splits=self.n_sample_splits
            ).split(X, self.y, groups=self.sample_array):
                train = np.intersect1d(sample_train, chr_train)
                test = np.intersect1d(sample_test, chr_test)

                # shuffle if requested
                if self.shuffle:
                    train = train.sample(frac=1, random_state=self.random_state)
                    test = test.sample(frac=1, random_state=self.random_state)

                yield train, test

    def get_n_splits(self):
        return self.n_splits


def count_loci(data: pd.DataFrame, label_col: str) -> int:
    "Count the number of loci in a dataset"
    return data.groupby("cell_id")[f"{label_col}_id"].nunique().sum()


class NpEncoder(json.JSONEncoder):
    "Encode pandas series and numpy arrays for json output"

    def default(self, obj):
        if isinstance(obj, pd.Series):
            return obj.to_numpy().tolist()
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


def parse_log(log: dict) -> pd.DataFrame:
    "Parse Model log file into a dataframe"
    out = []
    for f in log.keys():
        if isinstance(f, int):
            for s in ["train", "test"]:
                log[f][s].update(log[f]["model"])
                log[f][s]["fold"] = f
                out.append(log[f][s])

    return pd.DataFrame(out)


# TODO: get learning curves
class Model:
    def __init__(
        self,
        clf: AutoML,
        data: pd.DataFrame,
        features: list,
        label_col: str,
        rpm_filter: int,
        outfile: str,
    ) -> None:
        """
        Parameters
        ----------
        data : pd.DataFrame, data to fit model on
        features : list, feature columns to use for model
        label_col : str, column with labels
        rpm_filter : int, filter out windows with rpm < rpm_filter
        """

        # check inputs
        for c in [
            label_col,
            f"{label_col}_id",
            "cell_id",
            "Chromosome",
            "Start",
            "End",
            "rpm",
        ]:
            assert c in data.columns, f"data must have column {c}"
        for f in features:
            assert f in data.columns, f"data must have feature column {f}"
        assert data[label_col].dtype == bool, f"{label_col} must be boolean"

        # save attributes
        self.clf = clf
        self.total_pos_windows = data[label_col].sum()
        self.total_pos_loci = count_loci(data, label_col)
        self.data = data.loc[data["rpm"] >= rpm_filter, :].reset_index(drop=True)
        self.data["proba"] = np.nan  # initialize proba column
        self.filtered_pos_windows = self.data[label_col].sum()
        self.filtered_pos_loci = count_loci(self.data, label_col)
        self.features = features
        self.label_col = label_col
        self.rpm_filter = rpm_filter

        # check outfile
        assert not Path(
            outfile
        ).exists(), f"{outfile} already exists! Please choose a different name or delete the file if you want to overwrite it"
        self.outfile = outfile

        self.log = {
            "Input data": {
                "Total windows": len(data),
                "Total positive windows": self.total_pos_windows,
                "Total positive loci": self.total_pos_loci,
                "RPM threshold": rpm_filter,
                "Filtered windows": len(self.data),
                "Filtered positive windows": self.filtered_pos_windows,
                "Filtered positive loci": self.filtered_pos_loci,
                "Features": self.features,
                "Class label": self.label_col,
                "Number of cells": self.data.cell_id.nunique(),
                "Cells": self.data.cell_id.unique(),
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
        """

        # define train sets
        train = self.data.iloc[train_idx, :].reset_index()

        print(
            f"""
		Tuning model with {len(self.features)} features on {train.shape[0]} windows
		{train[self.label_col].sum()} positive windows ({count_loci(train, self.label_col)} loci)
		{len(train) - train[self.label_col].sum()} negative windows
		{train.Chromosome.nunique()} Chromosomes: {train.Chromosome.unique()}
		{train.donor_id.nunique()} Donors: {train.donor_id.unique()}
		{train.cell_id.nunique()} cells
		"""
        )

        # suppress UserWarning: `use_label_encoder` is deprecated in 1.7.0.
        warnings.filterwarnings("ignore", category=UserWarning)
        self.clf.fit(
            train[self.features],
            train[self.label_col],
            split_type=StratifiedGroupKFold(n_splits=n_hyper_splits),
            groups=train["Chromosome"],
            log_file_name=log_file_name,
        )

    def pred_eval(
        self, train_idx: np.ndarray, test_idx: np.ndarray
    ) -> dict[str, np.ndarray]:

        # predict and evaluate
        train = self.data.iloc[train_idx, :].reset_index()
        test = self.data.iloc[test_idx, :].reset_index()

        # make predictions for train and test sets
        result = []
        for name, df in {"train": train, "test": test}.items():
            proba = self.clf.predict_proba(df[self.features])[:, 1]

            # get the number of pos windows removed from this split due to rpm filter
            extra_pos_windows = (self.total_pos_windows - self.filtered_pos_windows) * (
                df[self.label_col].sum() / self.filtered_pos_windows
            )
            extra_pos_loci = (self.total_pos_loci - self.filtered_pos_loci) * (
                count_loci(df, self.label_col) / self.filtered_pos_loci
            )

            # get metrics
            metrics = Evaluator(df, proba, self.label_col).all_metrics(
                extra_pos_windows=extra_pos_windows, extra_pos_loci=extra_pos_loci
            )
            metrics["stage"] = name
            metrics["Chromosomes"] = df["Chromosome"].unique().tolist()
            metrics["Donors"] = df["donor_id"].unique().tolist()
            metrics["Cells"] = df["cell_id"].unique().tolist()
            metrics["n_windows"] = len(df)
            result.append(metrics.copy())

            if name == "test":
                # assign proba to self.data
                # convert test_idx to boolean
                # test_idx = np.zeros(len(self.data), dtype=bool)
                # test_idx[test.index] = True
                self.data.loc[test_idx, "proba"] = proba

        # add feature importances
        result.append(
            {
                "feature_importances": self.clf.feature_importances_,
                "features": self.clf.feature_names_in_,
                "best_estimator": self.clf.best_estimator,
                "best_config": self.clf.best_config,
                "total_loci_train": count_loci(train, self.label_col),
            }
        )

        return tuple(result)

    def cv(self, n_splits: int) -> None:
        # TODO: return predictions for each fold!

        # define splitter from cross-validation
        sgkf = StratifiedGroupKFold(n_splits=n_splits)

        # setup log file
        with open(self.outfile, "a") as f:
            for i, (train_idx, test_idx) in enumerate(
                sgkf.split(
                    self.data, self.data[self.label_col], groups=self.data["Chromosome"]
                )
            ):
                print(f"Fold {i+1}/{n_splits}")
                self.log[i + 1] = {}

                # tune model
                with tempfile.NamedTemporaryFile() as tmp:
                    self.tune(train_idx, n_hyper_splits=5, log_file_name=tmp.name)
                    (
                        time_history,
                        best_valid_loss_history,
                        _,
                        _,
                        _,
                    ) = get_output_from_log(
                        filename=tmp.name, time_budget=self.clf._settings["time_budget"]
                    )

                # predict and evaluate, saving results to log
                (
                    self.log[i + 1]["train"],
                    self.log[i + 1]["test"],
                    self.log[i + 1]["model"],
                ) = self.pred_eval(train_idx, test_idx)
                self.log[i + 1]["model"]["time_history"] = time_history
                self.log[i + 1]["model"][
                    "best_valid_loss_history"
                ] = best_valid_loss_history

        with open(self.outfile, "w") as f:
            json.dump(self.log, f, cls=NpEncoder)

    def get_results(self) -> pd.DataFrame:

        return parse_log(self.log)
