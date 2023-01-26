#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys, pickle, os
from gzip import GzipFile
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedGroupKFold
from imblearn.under_sampling import RandomUnderSampler

# get access to the src directory
sys.path.append((os.path.abspath("workflow")))
from src.model import Model

def read_prep_data(files):

    # read in data
    df = pd.concat([pd.read_parquet(fn) for fn in files]).set_index(
        [
            "chrom",
            "start",
            "end",
            "cell_id",
            "donor_id",
        ]
    )

    # make features
    features = [
        x
        for x in df.columns
        if not any(
            x.endswith(y)
            for y in [
                "chrom",
                "start",
                "end",
                "cell_id",
                "donor_id",
                "label",
                "build",
                "db",
            ]
        )
    ]

    # replace NA values with 0
    df[features] = df[features].fillna(0)

    # take minimum of features and 4e9 to avoid overflow error
    df[features] = np.minimum(df[features], 4e9)

    return df, features

if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    # read and prepare data
    df, features = read_prep_data(snakemake.input)

    # setup folds strategy
    kf = StratifiedGroupKFold(n_splits=snakemake.params.num_folds, shuffle = True, random_state=42)

    # setup test sampling strategy
    if snakemake.params.test_sampling_strategy:
        test_sampler = RandomUnderSampler(sampling_strategy=snakemake.params.test_sampling_strategy, random_state=42)
    else:
        test_sampler = None

    # create object for model
    m = Model(df, features, snakemake.params.min_reads)
    m = m.split(kf, snakemake.params.split_by, test_sampler)

    # save
    with GzipFile(snakemake.output[0], "wb") as f:
        pickle.dump(m, f)

    sys.stderr.close()
