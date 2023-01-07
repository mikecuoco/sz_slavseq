#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys
import pickle
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils import shuffle
from sklearn.metrics import precision_recall_curve, auc, confusion_matrix


def make_pipeline(clf_type, params):

    # TODO: add hyperparameter tuning
    if clf_type == "RandomForestClassifier":
        from sklearn.ensemble import RandomForestClassifier

        clf = RandomForestClassifier()
        pipe = Pipeline([(clf_type, clf)])
    elif clf_type == "LogisticRegression":
        from sklearn.linear_model import LogisticRegression

        clf = LogisticRegression()
        pipe = Pipeline([("scaler", StandardScaler()), (clf_type, clf)])
    elif clf_type == "SVC":
        from sklearn.svm import SVC

        clf = SVC()
        pipe = Pipeline([("scaler", StandardScaler()), (clf_type, clf)])
    elif clf_type == "MLPClassifier":
        from sklearn.neural_network import MLPClassifier

        clf = MLPClassifier()
        pipe = Pipeline([(clf_type, clf)])
    elif clf_type == "CascadeForestClassifier":
        from deepforest import CascadeForestClassifier

        clf = CascadeForestClassifier()
        pipe = Pipeline([(clf_type, clf)])
    elif clf_type == "DecisionTreeClassifier":
        from sklearn.tree import DecisionTreeClassifier

        clf = DecisionTreeClassifier()
        pipe = Pipeline([(clf_type, clf)])
    elif clf_type == "GradientBoostingClassifier":
        from sklearn.ensemble import GradientBoostingClassifier

        clf = GradientBoostingClassifier()
        pipe = Pipeline([(clf_type, clf)])
    elif clf_type == "HistGradientBoostingClassifier":
        from sklearn.ensemble import HistGradientBoostingClassifier

        clf = HistGradientBoostingClassifier()
        pipe = Pipeline([(clf_type, clf)])
    elif clf_type == "XGBClassifier":
        from xgboost import XGBClassifier

        clf = XGBClassifier()
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
        y_true = pickle.load(f)
    with open(snakemake.input.label_encoder, "rb") as f:
        le = pickle.load(f)

    assert set(y_true.keys()) == set(
        features.keys()
    ), "features and labels must have the same number of folds"

    # train/test/metrics for each fold in a loop
    y_pred, y_proba, metrics = {}, {}, {}
    for fold in features.keys():
        y_pred[fold], y_proba[fold], metrics[fold] = {}, {}, {}

        # train
        pipe = make_pipeline(snakemake.params.model_name, snakemake.params.model_params)
        pipe.fit(features[fold]["train"], y_true[fold]["train"])

        # make shuffled test data
        y_true[fold]["test_shuffled"] = shuffle(y_true[fold]["test"], random_state=42)

        for stage in ["train", "test", "test_shuffled"]:
            metrics[fold][stage] = {}

            # make predictions
            if stage == "test_shuffled":
                y_pred[fold][stage] = pipe.predict(features[fold]["test"])
                y_proba[fold][stage] = pipe.predict_proba(features[fold]["test"])
            else:
                y_pred[fold][stage] = pipe.predict(features[fold][stage])
                y_proba[fold][stage] = pipe.predict_proba(features[fold][stage])

            # get confusion matrices
            # raw
            metrics[fold][stage]["cm"] = confusion_matrix(
                y_true=le.inverse_transform(y_true[fold][stage]),
                y_pred=le.inverse_transform(y_pred[fold][stage]),
                labels=le.classes_,
            )

            # normalized
            metrics[fold][stage]["cm_norm"] = confusion_matrix(
                y_true=le.inverse_transform(y_true[fold][stage]),
                y_pred=le.inverse_transform(y_pred[fold][stage]),
                labels=le.classes_,
                normalize="true",
            )

            for label in ["KNRGL", "RL1", "OTHER"]:
                metrics[fold][stage][label] = {}

                # get prcurve
                pr = pd.DataFrame()
                pr["precision"], pr["recall"], _ = precision_recall_curve(
                    y_true[fold][stage],
                    y_proba[fold][stage][:, le.transform([label])[0]],
                    pos_label=le.transform([label])[0],
                )
                metrics[fold][stage][label]["prcurve"] = pr
                metrics[fold][stage][label]["auprc"] = auc(
                    pr["recall"], pr["precision"]
                )

    # save predictions
    with open(snakemake.output.proba, "wb") as f:
        pickle.dump(y_proba, f)
    with open(snakemake.output.pred, "wb") as f:
        pickle.dump(y_pred, f)

    # save metrics
    with open(snakemake.output.metrics, "wb") as f:
        pickle.dump(metrics, f)

    sys.stderr.close()
