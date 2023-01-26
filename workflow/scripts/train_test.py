#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys, pickle, os
from gzip import GzipFile
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from imblearn.under_sampling import RandomUnderSampler
# get access to the src directory
sys.path.append((os.path.abspath("workflow")))

def make_pipeline(clf_type, params):

    # TODO: add hyperparameter tuning
    if clf_type == "KNeighborsClassifier":
        from sklearn.neighbors import KNeighborsClassifier

        clf = KNeighborsClassifier()
        pipe = [("scaler", StandardScaler())]

    elif clf_type == "DecisionTreeClassifier":
        from sklearn.tree import DecisionTreeClassifier

        clf = DecisionTreeClassifier()

    elif clf_type == "RandomForestClassifier":
        from sklearn.ensemble import RandomForestClassifier

        clf = RandomForestClassifier()

    elif clf_type == "BalancedRandomForestClassifier":
        from imblearn.ensemble import BalancedRandomForestClassifier

        clf = BalancedRandomForestClassifier()

    elif clf_type == "CascadeForestClassifier":
        from deepforest import CascadeForestClassifier

        clf = CascadeForestClassifier()

    elif clf_type == "GradientBoostingClassifier":
        from sklearn.ensemble import GradientBoostingClassifier

        clf = GradientBoostingClassifier()

    elif clf_type == "HistGradientBoostingClassifier":
        from sklearn.ensemble import HistGradientBoostingClassifier

        clf = HistGradientBoostingClassifier()

    elif clf_type == "AdaBoostClassifier":
        from sklearn.ensemble import AdaBoostClassifier

        clf = AdaBoostClassifier()

    elif clf_type == "XGBClassifier":
        from xgboost import XGBClassifier

        clf = XGBClassifier()

    elif clf_type == "LogisticRegression":
        from sklearn.linear_model import LogisticRegression

        clf = LogisticRegression()
        pipe = [("scaler", StandardScaler())]

    elif clf_type == "SVC":
        from sklearn.svm import SVC

        clf = SVC()
        pipe = [("scaler", StandardScaler())]

    elif clf_type == "MLPClassifier":
        from sklearn.neural_network import MLPClassifier

        clf = MLPClassifier()

    if "pipe" in locals():
        pipe.append((clf_type, clf))
        pipe = Pipeline(pipe)
    else:
        pipe = Pipeline([(clf_type, clf)])

    # add user-specified hyperparameters
    if params != None:
        params = {f"{clf_type}__{k}": v for k, v in params.items()}
        pipe.set_params(**params)

    return pipe


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    with GzipFile(snakemake.input[0], "rb") as f:
        m = pickle.load(f)

    # setup
    if snakemake.params.train_sampling_strategy:
        train_sampler = RandomUnderSampler(sampling_strategy=snakemake.params.train_sampling_strategy, random_state=42)
    else:
        train_sampler = None

    clf = make_pipeline(snakemake.params.model_name, snakemake.params.model_params)
    m = m.train_test(clf, train_sampler, snakemake.wildcards.model_id)

    with GzipFile(snakemake.output[0], "wb") as f:
        pickle.dump(m, f)

    sys.stderr.close()
