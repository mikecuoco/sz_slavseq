#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys
import pandas as pd
import pickle
from sklearn.utils import shuffle
from sklearn.metrics import precision_recall_curve


sys.stderr = open(snakemake.log[0], "w")

# get label encoder
with open(snakemake.input.label_encoder, "rb") as f:
    le = pickle.load(f)

# get model results
with open(snakemake.input.proba, "rb") as f:
    proba = pickle.load(f)

# get labels
with open(snakemake.input.labels, "rb") as f:
    y = pickle.load(f)

# get PR curves in a loop
df_list = []
for fold in proba.keys():
    for stage in ["train", "test"]:
        for label in ["KNRGL", "RL1", "OTHER"]:

            # TODO: figure out how to return thresholds
            _df = pd.DataFrame()
            _df["precision"], _df["recall"], _ = precision_recall_curve(
                y[fold][stage],
                proba[fold][stage][:, le.transform([label])[0]].round(2),
                pos_label=le.transform([label])[0],
            )
            _df["stage"], _df["label"], _df["fold"] = stage, label, fold + 1
            _df["model_id"] = snakemake.wildcards.model_id
            df_list.append(_df)

            if stage == "test":
                _df = pd.DataFrame()
                _df["precision"], _df["recall"], _ = precision_recall_curve(
                    shuffle(y[fold][stage]),
                    proba[fold][stage][:, le.transform([label])[0]].round(2),
                    pos_label=le.transform([label])[0],
                )
                _df["stage"], _df["label"], _df["fold"] = (
                    "test_shuffled",
                    label,
                    fold + 1,
                )
                df_list.append(_df)


df = pd.concat(df_list).reset_index(drop=True)
df.to_pickle(snakemake.output[0])

sys.stderr.close()
