#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = 'Michael Cuoco'

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.stderr = open(snakemake.log[0], "w")

df = pd.concat([pd.read_pickle(f) for f in snakemake.input]).round(1)

# Plot PR curves
sns.set_style("ticks")
fig = sns.relplot(
    data=df,
    x="recall",
    y="precision",
    hue="label",
    col="stage",
    row="model_id",
    kind="line",
    markers=True,
)

plt.savefig(snakemake.output[0], format="png", dpi=200)
plt.clf()

# TODO: compute confusion matrix
# TODO: report feature importance