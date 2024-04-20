#!/usr/bin/env python
# Created on: Dec 14, 2023 at 4:37:45 PM
__author__ = "Michael Cuoco"


import pickle
import pandas as pd

# ML
from flaml import AutoML
from flaml import tune
from sklearn.model_selection import StratifiedGroupKFold

# logging
import logging
from flaml.automl.logger import logger
from flaml.tune.tune import logger as tune_logger


# define custom hp
custom_hp = {
    "xgboost": {
        "n_estimators": {
            "domain": tune.randint(lower=100, upper=500),
            "low_cost_init_value": 250,
        },
        "max_leaves": {
            "domain": tune.lograndint(lower=6, upper=1024),
            "low_cost_init_value": 50,
        },
        "max_depth": {
            "domain": tune.randint(lower=3, upper=100),
            "low_cost_init_value": 5,
        },
    }
}

# read data
if __name__ == "__main__":
    fh = logging.FileHandler(snakemake.log[0], mode="w")  # type: ignore
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    tune_logger.addHandler(fh)

    data = pd.read_parquet(snakemake.input.data)  # type: ignore

    # define features
    features = snakemake.config["features"]  # type: ignore

    for f in features:
        assert f in data.columns, f"Feature {f} not in dataframe columns"
    assert len(features) == len(set(features)), "Duplicate features"

    # define model
    automl = AutoML(
        task="binary",
        estimator_list=["xgboost"],
        metric="ap",
        eval_method="cv",
        max_iter=snakemake.params.max_iter,  # type: ignore
        n_jobs=snakemake.threads,  # type: ignore
        skip_transform=False,  # don't preprocess data
        auto_augment=False,  # don't augment rare classes
        early_stop=True,
        retrain_full=False,
        verbose=4,
        seed=123,
        log_training_metric=True,
        custom_hp=custom_hp,
        log_file_name=snakemake.output.history,  # type: ignore
    )

    # tune model
    sgkf = StratifiedGroupKFold(
        n_splits=5,
        shuffle=True,
        random_state=snakemake.params.random_state,  # type: ignore
    )

    automl.fit(
        X_train=data[features],
        y_train=data["KNRGL"],
        groups=data["Chromosome"],
        split_type=sgkf,
    )

    # save model and hyperparameters
    clf = automl.model.estimator
    with open(snakemake.output.model, "wb") as f:  # type: ignore
        pickle.dump(clf, f, pickle.HIGHEST_PROTOCOL)
    automl.save_best_config(snakemake.output.best_hp)  # type: ignore
