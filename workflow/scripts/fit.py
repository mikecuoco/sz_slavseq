#!/usr/bin/env python
# Created on: May 1, 2023 at 5:04:13 PM
__author__ = "Michael Cuoco"

import sys

# set the log file
sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import pickle as pkl
from flaml import AutoML
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold, GridSearchCV, cross_val_score
from sklearn.metrics import make_scorer
import pdb

# read the data
df = pd.concat([pd.read_parquet(f) for f in snakemake.input])
assert set(df["label"].unique()) == {"KNRGL", "OTHER"}, "labels are not correct"

# save all features for model to list
features = [
    x
    for x in df.columns
    if not any(
        x.endswith(y)
        for y in [
            "Chromosome",
            "Start",
            "End",
            "cell_id",
            "donor_id",
            "label",
            "blacklist",
            "index",
            "reads",
            "fwd",
            "rev",
            "reads_bg",
            "fwd_bg",
            "rev_bg",
        ]
    )
]

assert all(x in df.columns for x in features), "some features are not in the dataframe"

# # encode labels
print("encoding labels as integers...")
df["label_encoded"] = df["label"].map({"KNRGL": 1, "OTHER": 0})


# fit the model
clf = AutoML()
clf.fit(
    X_train=df[features],
    y_train=df["label_encoded"],
    task="classification",
    metric="f1",
    n_jobs=snakemake.threads,
    estimator_list=["xgboost"],
    groups=df["cell_id"],
    early_stop=True,
    time_budget=600,
    log_file_name=sys.stderr.name,
)

# save the model
with open(snakemake.output["model"], "wb") as f:
    pkl.dump(clf, f, pkl.HIGHEST_PROTOCOL)

# save the features
with open(snakemake.output["features"], "w") as f:
    f.write("\n".join(features))

sys.stderr.close()
