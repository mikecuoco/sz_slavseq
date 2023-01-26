from time import time
import sys
import numpy as np
import pandas as pd


class Model:
    def __init__(self, clf: object, train_sampler=None, id: str = None):
        self.clf_ = clf
        self.train_sampler_ = train_sampler
        self.id_ = id

    def train_test(self, folds: dict, features: list):

        self.res_ = {}
        self.features_ = features

        for f in folds.keys():
            self.res_[f] = {}

            # resample the training set if specified
            if self.train_sampler_:
                folds[f]["train"], _ = self.train_sampler_.fit_resample(
                    folds[f]["train"].loc[
                        :,
                    ],
                    folds[f]["train"].loc[:, "label"],
                )

            # fit the model
            self.res_[f]["clf"] = self.clf_.fit(
                folds[f]["train"].loc[:, features],
                folds[f]["train"].loc[:, "label"],
            )

            for stage in ["train", "test"]:
                self.res_[f][stage] = {}
                self.res_[f][stage]["label"] = folds[f][stage].loc[:, "label"]

                self.res_[f][stage]["pred"] = self.res_[f]["clf"].predict(
                    folds[f][stage].loc[:, features]
                )

                proba = self.res_[f]["clf"].predict_proba(
                    folds[f][stage].loc[:, features]
                )

                for label in self.res_[f]["clf"].classes_:
                    (i,) = np.where(
                        self.res_[f]["clf"].classes_ == label
                    )  # get index of label
                    self.res_[f][stage][f"proba_{label}"] = proba[:, i]

        return self
