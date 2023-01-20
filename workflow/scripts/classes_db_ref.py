#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# TODO: color by donor race/schizophrenia status
donors = pd.concat([pd.read_pickle(fn).reset_index() for fn in snakemake.input])

plot_df = (
    donors.value_counts(["label", "build", "db", "cell_id", "donor_id"])
    .to_frame("count")
    .reset_index()
)
plot_df["ref_db"] = plot_df["build"] + "_" + plot_df["db"].astype(str)
plot_df = plot_df[plot_df["label"] == "KNRGL"]

sns.boxplot(plot_df, x="count", y="ref_db", hue="donor_id")
plt.savefig(snakemake.output[0], format="png", dpi=200, bbox_inches="tight")
