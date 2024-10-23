#!/usr/bin/env python
# Created on: Nov 8, 2023 at 9:39:33 PM
__author__ = "Michael Cuoco"


from flaml.automl.model import (
    SKLearnEstimator,
    XGBoostEstimator,
    LGBMEstimator,
    XGBoostSklearnEstimator,
)
from flaml import tune
from sklearn import linear_model, pipeline
from sklearn.preprocessing import StandardScaler


class ElasticNetLR(SKLearnEstimator):
    def __init__(self, task="classification", n_jobs=None, **config):
        super().__init__(task, **config)

        self.estimator_class = self._create_pipeline

    @staticmethod
    def _create_pipeline(**kwargs):
        lr = linear_model.LogisticRegression(penalty="elasticnet", **kwargs)
        return pipeline.Pipeline([("scaler", StandardScaler()), ("lr", lr)])

    @classmethod
    def search_space(cls, data_size, task):
        space = {
            "C": {
                "domain": tune.loguniform(lower=1, upper=10),
                "init_value": 1,
            },
            "l1_ratio": {
                "domain": tune.loguniform(lower=0, upper=1),
                "init_value": 0.5,
            },
        }
        return space


def MyXGB(XGBoostEstimator):
    @classmethod
    def search_space(cls, data_size, task):
        upper = min(32768, int(data_size))
        return {
            "n_estimators": {
                "domain": tune.randint(lower=50, upper=5000),
                "low_cost_init_value": 50,
            },
            "max_leaves": {
                "domain": tune.lograndint(lower=6, upper=512),
                "low_cost_init_value": 50,
            },
        }


def MyLGBM(LGBMEstimator):
    @classmethod
    def search_space(cls, data_size, task):
        upper = min(32768, int(data_size))
        return {
            "n_estimators": {
                "domain": tune.randint(lower=50, upper=5000),
                "low_cost_init_value": 50,
            },
            "num_leaves": {
                "domain": tune.lograndint(lower=6, upper=512),
                "low_cost_init_value": 50,
            },
        }
