#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys
import pickle
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from deepforest import CascadeForestClassifier


def make_pipeline(clf_type, params):

    # TODO: add hyperparameter tuning
    # TODO: add DecisionTreeClassifier
    # TODO: add XGBoost
    # TODO: add LightGBM
    # TODO: add GradientBoostingClassifier
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

    with open(snakemake.input.features, "rb") as f:
        features = pickle.load(f)
    with open(snakemake.input.labels, "rb") as f:
        labels = pickle.load(f)

    assert set(labels.keys()) == set(
        features.keys()
    ), "features and labels must have the same number of folds"

    pred, proba = {}, {}

    for fold in features.keys():
        pred[fold], proba[fold] = {}, {}

        # train/test
        pipe = make_pipeline(snakemake.params.model_name, snakemake.params.model_params)
        pipe.fit(features[fold]["train"], labels[fold]["train"])
        pred[fold]["train"] = pipe.predict(features[fold]["train"])
        pred[fold]["test"] = pipe.predict(features[fold]["test"])
        proba[fold]["train"] = pipe.predict_proba(features[fold]["train"])
        proba[fold]["test"] = pipe.predict_proba(features[fold]["test"])

    # save predictions
    with open(snakemake.output.proba, "wb") as f:
        pickle.dump(proba, f)
    with open(snakemake.output.pred, "wb") as f:
        pickle.dump(pred, f)

    sys.stderr.close()
