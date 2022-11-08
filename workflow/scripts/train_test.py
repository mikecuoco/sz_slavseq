#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys
import pickle
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from deepforest import CascadeForestClassifier


def make_pipeline(clf_type, params):

    # TODO: add hyperparameter tuning
    if clf_type == "RandomForestClassifier":
        clf = RandomForestClassifier()
        pipe = Pipeline([(clf_type, clf)])
    elif clf_type == "LogisticRegression":
        clf = LogisticRegression()
        pipe = Pipeline([("scaler", StandardScaler()), (clf_type, clf)])
    elif clf_type == "SVC":
        clf = SVC()
        pipe = Pipeline([("scaler", StandardScaler()), (clf_type, clf)])
    elif clf_type == "MLPClassifier":
        clf = MLPClassifier()
        pipe = Pipeline([(clf_type, clf)])
    elif clf_type == "CascadeForestClassifier":
        clf = CascadeForestClassifier()
        pipe = Pipeline([(clf_type, clf)])

    # add user-specified hyperparameters
    if params != None:
        params = {f"{clf_type}__{k}": v for k, v in params.items()}
        pipe.set_params(**params)

    return pipe


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    for fold in range(0, snakemake.params.num_folds):
        pipe = make_pipeline(snakemake.params.model_name, snakemake.params.model_params)

        # train
        X_train = pd.read_pickle(snakemake.input.train_features[fold])
        y_train = pd.read_pickle(snakemake.input.train_labels[fold])
        pipe.fit(X_train, y_train)

        y_train_pred = pipe.predict(X_train)
        with open(snakemake.output.train_pred[fold], "wb") as f:
            pickle.dump(y_train_pred, f)
        y_train_proba = pipe.predict_proba(X_train)
        with open(snakemake.output.train_proba[fold], "wb") as f:
            pickle.dump(y_train_proba, f)

        # test
        X_test = pd.read_pickle(snakemake.input.test_features[fold])
        y_test_pred = pipe.predict(X_test)
        with open(snakemake.output.test_pred[fold], "wb") as f:
            pickle.dump(y_test_pred, f)
        y_test_proba = pipe.predict_proba(X_test)
        with open(snakemake.output.test_proba[fold], "wb") as f:
            pickle.dump(y_test_proba, f)

    # TODO: train model on all data and return that model

    sys.stderr.close()
