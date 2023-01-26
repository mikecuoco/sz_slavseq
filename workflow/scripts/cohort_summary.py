#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import pandas as pd
import seaborn as sns

# read in data
df = pd.read_csv(snakemake.input[0], sep="\t", dtype={"sample": str, "donor": str})

# donors
plot_df = df[["donor", "age", "race", "diagnosis"]].drop_duplicates()
fig = sns.swarmplot(plot_df, y="age", x="diagnosis", hue="race")

# cells per tissue/donor
plot_df = df.value_counts(["donor", "region"]).to_frame("cells").reset_index()
fig = sns.barplot(data=plot_df, y="cells", x="donor", hue="region")

# TODO: reads per cell per donor per dna_type
