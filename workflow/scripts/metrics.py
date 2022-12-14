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
df_list = []
for fold in range(0, snakemake.params.num_folds):
    for stage in ["train", "test"]:
        with open(snakemake.input[f"{stage}_labels"][fold], "rb") as f:
            y = pickle.load(f)

        file = [f for f in snakemake.input[f"{stage}_proba"]][fold]
        with open(file, "rb") as f:
            y_proba = pickle.load(f)

        for label in ["KNRGL","RL1","OTHER"]:
            _df = pd.DataFrame()
            
            _df["precision"], _df["recall"], _ = precision_recall_curve(
                y, y_proba[:, le.transform([label])[0]], pos_label=le.transform([label])[0]
            )
            
            _df["stage"] = stage
            _df["label"] = label
            _df["fold"] = fold + 1
            df_list.append(_df)
        
df = pd.concat(df_list).reset_index(drop=True)

# Plot PR curves
sns.set_style("ticks")
fig = sns.relplot(
    data=df,
    x="recall",
    y="precision",
    hue="label",
    col="stage",
    row="fold",
    kind="line",
    markers=True,
)
fig.set(
    ylim=[0, 1],
    xlim=[0, 1]
)
plt.savefig(snakemake.output.prcurve, format="png", dpi=300)

# TODO: compute confusion matrix
# TODO: report feature importance

sys.stderr.close()
