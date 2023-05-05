#!/usr/bin/env python
# Created on: May 1, 2023 at 4:53:24 PM
__author__ = "Michael Cuoco"

import sys

# set the log file
sys.stderr = open(snakemake.log[0], "w")

import pickle as pkl
import pandas as pd
import pyranges as pr
from get_labels import filter_windows

# load the model
with open(snakemake.input["model"], "rb") as f:
    clf = pkl.load(f)

# load the data
df = pd.read_parquet(snakemake.input["data"][0])

# load the blacklist and segdups
sv_blacklist = pr.read_bed(snakemake.input.sv_blacklist[0]).df
segdups = pr.read_bed(snakemake.input.segdups[0]).df

# filter the data based on dictionary of boolean filters
print("filtering data...", file=sys.stderr)
df = filter_windows(df, sv_blacklist, segdups)
assert len(df) > 0, "no windows left after filtering"

# make predictions
print("making predictions...", file=sys.stderr)
pred = clf.predict(df[clf.feature_names_in_])
df["pred"] = ["KNRGL" if x == 1 else "OTHER" for x in pred]
df["KNRGL_proba"] = clf.predict_proba(df[clf.feature_names_in_])[:, 1]

# save the predictions
print("saving predictions...", file=sys.stderr)
df.to_parquet(snakemake.output[0])

sys.stderr.close()
