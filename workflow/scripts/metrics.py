#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys
import pandas as pd
import pickle
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
import seaborn as sns

sys.stderr = open(snakemake.log[0], "w")

# get label encoder
with open(snakemake.input.label_encoder, "rb") as f:
    le = pickle.load(f)

# get model results in loop
df = pd.DataFrame()
for fold in range(0, snakemake.params.num_folds):
    for stage in ["train", "test"]:
        file = snakemake.input[f"{stage}_labels"][fold]
        with open(file, "rb") as f:
            y = pickle.load(f)

        for model in snakemake.params.models:
            _df = pd.DataFrame()

            file = [f for f in snakemake.input[f"{stage}_proba"] if model in f][fold]

            with open(file, "rb") as f:
                y_proba = pickle.load(f)[:, le.transform(['KNRGL'])[0]]

            _df["precision"], _df["recall"], _ = precision_recall_curve(
                y, y_proba, pos_label=le.transform(["KNRGL"])[0]
            )

            _df["stage"] = stage
            _df["model"] = model
            _df["fold"] = fold + 1
            df = pd.concat([df, _df]).reset_index(drop=True)

# Plot PR curves
sns.set_style("ticks")
fig = sns.relplot(
    data=df,
    x="recall",
    y="precision",
    hue="model",
    col="stage",
    kind="line",
    markers=True,
)
fig.set(
    ylim=[0, 1],
    xlim=[0, 1],
    title="Precision-Recall Curve: KNRGL vs. RL1, OTHER",
)
plt.savefig(snakemake.output.prcurve, format="svg", dpi=300)

# TODO: compute confusion matrix
# TODO: report feature importance

sys.stderr.close()
