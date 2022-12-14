#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = 'Michael Cuoco'

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.concat([pd.read_pickle(f) for f in snakemake.input]).reset_index()

# classes per cell
# TODO: move legend to outside of plot
plot_df = df.value_counts(["label","donor_id","cell_id"]).to_frame('count').reset_index()
fig = sns.boxplot(plot_df, x="count", y="donor_id", hue="label")
fig.set_xscale("log")
plt.savefig(snakemake.output.classes_per_cell, format="png", dpi=300)

# classes per donor
# TODO: move legend to outside of plot
plot_df = df.value_counts(["label","donor_id"]).to_frame('count').reset_index()
fig = sns.barplot(plot_df, x="count", y="donor_id", hue="label")
fig.set_xscale("log")
plt.savefig(snakemake.output.classes_per_donor, format="png", dpi=300)

# features per class
plot_df = df.melt(id_vars=["chrom","start","end","donor_id","cell_id","label"], var_name="feature")
fig = sns.FacetGrid(plot_df, col="feature", col_wrap=3, sharex=False)
fig.map_dataframe(sns.boxplot, x="value", y="label", hue="label")
plt.savefig(snakemake.output.features_per_class, format="png", dpi=300)

# TODO: evaluate folds