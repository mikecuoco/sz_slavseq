#!/usr/bin/env python
# Created on: May 1, 2023 at 5:04:13 PM
__author__ = "Michael Cuoco"

import sys
from tempfile import NamedTemporaryFile
import pyarrow.parquet as pq
import pandas as pd
import pickle as pkl
from flaml import AutoML
from flaml.automl.data import get_output_from_log
from sklearn import model_selection
from sklearn.metrics import (
    precision_recall_curve,
    roc_curve,
    PrecisionRecallDisplay,
    RocCurveDisplay,
)
import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":
    # set the log file
    sys.stderr = open(snakemake.log[0], "w")

    # read the data
    df = pd.concat([pq.read_table(f).to_pandas() for f in snakemake.input])

    # encode labels as integers
    df["label_encoded"] = df["label"].map({"KNRGL": 1, "OTHER": 0})

    # define features
    features = []
    keys = ["_q", "frac", "gini", "bias"]
    for c in df.columns:
        for k in keys:
            if k in c:
                features.append(c)

    # define model training settings
    flaml_settings = dict(
        task="classification",
        n_jobs=snakemake.threads,
        estimator_list=["xgboost", "rf"],
        early_stop=True,
        skip_transform=True,  # don't preprocess data
        auto_augment=False,  # don't augment rare classes
        time_budget=600,
    )

    m = Model(df, features, val_chr="chr1")
    m.fit(df, features, **flaml_settings)

    # save the model
    with open(snakemake.output["model"], "wb") as f:
        pkl.dump(clf, f, pkl.HIGHEST_PROTOCOL)

    # save the features
    with open(snakemake.output["features"], "w") as f:
        f.write("\n".join(features))

    sys.stderr.close()
