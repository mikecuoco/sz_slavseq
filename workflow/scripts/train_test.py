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
import pdb

def make_pipeline(clf_type):

    if clf_type == "RandomForestClassifier":
        clf = RandomForestClassifier(n_estimators=100, bootstrap=True, oob_score=True)
        pipe = Pipeline([(clf_type, clf)])
    elif clf_type == "LogisticRegression":
        clf = LogisticRegression(solver='saga', penalty='elasticnet', l1_ratio=0.5, max_iter=1e6)
        pipe = Pipeline([("scaler", StandardScaler()), (clf_type, clf)])
    elif clf_type == "SVC":
        clf = SVC(probability=True)
        pipe = Pipeline([("scaler", StandardScaler()), (clf_type, clf)])
    elif clf_type == "MLPClassifier":
        clf = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(5, 2), random_state=1, max_iter=1e6)
        pipe = Pipeline([(clf_type, clf)])

    return pipe


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    pipe = make_pipeline(snakemake.wildcards.model)

    for fold in range(0, snakemake.params.num_folds):
        # train
        X_train = pd.read_pickle(snakemake.input.train_features[fold])
        y_train = pd.read_pickle(snakemake.input.train_labels[fold])
        pipe.fit(X_train, y_train)

        y_train_pred = pipe.predict(X_train)
        pickle.dump(y_train_pred, open(snakemake.output.train_predictions[fold], "wb"))
        y_train_proba = pipe.predict_proba(X_train)
        pickle.dump(y_train_proba, open(snakemake.output.train_probabilities[fold], "wb"))

        # test
        X_test = pd.read_pickle(snakemake.input.test_features[fold])
        y_test_pred = pipe.predict(X_test)
        pickle.dump(y_test_pred, open(snakemake.output.test_predictions[fold], "wb"))
        y_test_proba = pipe.predict_proba(X_test)
        pickle.dump(y_test_proba, open(snakemake.output.test_probabilities[fold], "wb"))

    sys.stderr.close()