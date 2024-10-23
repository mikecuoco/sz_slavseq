"""Custom FLAML estimators."""

from flaml.automl.model import SKLearnEstimator
from flaml.automl.model import LGBMEstimator
from flaml import tune
from sklearn import linear_model
from sklearn import pipeline
from sklearn.preprocessing import StandardScaler

from partitioned_ensemble import PartitionedEnsembleRows, PartitionedEnsembleColumns


class Lasso(SKLearnEstimator):
    def __init__(self, task="regression", n_jobs=None, **config):
        super().__init__(task, **config)

        if task not in ["regression"]:
            raise ValueError("Lasso model is only for regression tasks")

        self.estimator_class = linear_model.Lasso

    @classmethod
    def search_space(cls, data_size, task):
        space = {
            "alpha": {
                "domain": tune.loguniform(lower=1e-10, upper=2.0),
                "init_value": 1e-4,
            },
            "max_iter": {
                "domain": tune.lograndint(lower=800, upper=10000),
                "init_value": 1000,
                "low_cost_init_value": 1000,
            },
            "tol": {
                "domain": tune.loguniform(lower=1e-7, upper=5e-3),
                "init_value": 1e-4,
            },
            "selection": {
                "domain": tune.choice(["cyclic", "random"]),
                "init_value": "cyclic",
            },
        }
        return space


class RowPartitionedLasso(SKLearnEstimator):
    """Lasso where data is first partitioned into n parts by sample and
    a model is trained on each part.
    """

    def __init__(self, task="regression", n_jobs=None, **config):
        super().__init__(task, estimator=linear_model.Lasso(), n_partitions=3, **config)

        if task not in ["regression"]:
            raise ValueError("RowPartitionedLasso model is only for regression tasks")

        self.estimator_class = PartitionedEnsembleRows

    @classmethod
    def search_space(cls, data_size, task):
        space = {
            "alpha": {
                "domain": tune.loguniform(lower=1e-10, upper=2.0),
                "init_value": 1e-4,
            },
            "max_iter": {
                "domain": tune.lograndint(lower=800, upper=10000),
                "init_value": 1000,
                "low_cost_init_value": 1000,
            },
            "tol": {
                "domain": tune.loguniform(lower=1e-7, upper=5e-3),
                "init_value": 1e-4,
            },
            "selection": {
                "domain": tune.choice(["cyclic", "random"]),
                "init_value": "cyclic",
            },
        }
        return space


class ElasticNet(SKLearnEstimator):
    def __init__(self, task="regression", n_jobs=None, **config):
        super().__init__(task, **config)

        if task not in ["regression"]:
            raise ValueError("ElasticNet model is only for regression tasks")

        self.estimator_class = linear_model.ElasticNet

    @classmethod
    def search_space(cls, data_size, task):
        space = {
            "alpha": {
                "domain": tune.loguniform(lower=1e-10, upper=2.0),
                "init_value": 1e-4,
            },
            "l1_ratio": {
                "domain": tune.uniform(0.01, 1),
                "init_value": 0.5,
            },
            "max_iter": {
                "domain": tune.lograndint(lower=800, upper=10000),
                "init_value": 1000,
                "low_cost_init_value": 1000,
            },
            "tol": {
                "domain": tune.loguniform(lower=1e-7, upper=5e-3),
                "init_value": 1e-4,
            },
            "selection": {
                "domain": tune.choice(["cyclic", "random"]),
                "init_value": "cyclic",
            },
        }
        return space


class RowPartitionedElasticNet(SKLearnEstimator):
    """Lasso where data is first partitioned into n parts by sample and
    a model is trained on each part.
    """

    def __init__(self, task="regression", n_jobs=None, **config):
        super().__init__(
            task, estimator=linear_model.ElasticNet(), n_partitions=3, **config
        )

        if task not in ["regression"]:
            raise ValueError("RowPartitionedLasso model is only for regression tasks")

        self.estimator_class = PartitionedEnsembleRows

    @classmethod
    def search_space(cls, data_size, task):
        space = {
            "alpha": {
                "domain": tune.loguniform(lower=1e-10, upper=2.0),
                "init_value": 1e-4,
            },
            "l1_ratio": {
                "domain": tune.uniform(0.01, 1),
                "init_value": 0.5,
            },
            "max_iter": {
                "domain": tune.lograndint(lower=800, upper=10000),
                "init_value": 1000,
                "low_cost_init_value": 1000,
            },
            "tol": {
                "domain": tune.loguniform(lower=1e-7, upper=5e-3),
                "init_value": 1e-4,
            },
            "selection": {
                "domain": tune.choice(["cyclic", "random"]),
                "init_value": "cyclic",
            },
        }
        return space


class OrthogonalMatchingPursuit(SKLearnEstimator):
    def __init__(self, task="regression", n_jobs=None, sample_weight=None, **config):
        super().__init__(task, normalize=False, **config)

        if task not in ["regression"]:
            raise ValueError(
                "OrthogonalMatchingPursuit model is only for regression tasks"
            )

        self.estimator_class = linear_model.OrthogonalMatchingPursuit

    @classmethod
    def search_space(cls, data_size, task):
        n_variants = data_size[1]
        space = {
            "n_nonzero_coefs": {
                "domain": tune.lograndint(1, upper=n_variants * 0.01),
                "init_value": max(int(n_variants * 0.001), 1),
            }
        }
        return space


class OrthogonalMatchingPursuitNormalized(SKLearnEstimator):
    def __init__(self, task="regression", n_jobs=None, sample_weight=None, **config):
        super().__init__(task, normalize=False, **config)

        if task not in ["regression"]:
            raise ValueError(
                "OrthogonalMatchingPursuit model is only for regression tasks"
            )

        # Create a pipeline that includes the StandardScaler and
        # OrthogonalMatchingPursuit
        self.estimator_class = self._create_pipeline

    @staticmethod
    def _create_pipeline(**kwargs):
        omp = linear_model.OrthogonalMatchingPursuit(**kwargs)
        return pipeline.Pipeline([("scaler", StandardScaler()), ("omp", omp)])

    @classmethod
    def search_space(cls, data_size, task):
        n_variants = data_size[1]
        space = {
            "n_nonzero_coefs": {
                "domain": tune.lograndint(1, upper=n_variants * 0.01),
                "init_value": max(int(n_variants * 0.001), 1),
            }
        }
        return space


class LGBMEstimatorHigherFloor(LGBMEstimator):
    """LGBMEstimator with higher floor for n_estimators."""

    @classmethod
    def search_space(cls, data_size, **params):
        return {
            "n_estimators": {
                "domain": tune.randint(lower=500, upper=5000),
                "init_value": 1000,
            },
            "num_leaves": {
                "domain": tune.lograndint(lower=6, upper=512),
                "init_value": 31,
            },
            "min_child_samples": {
                "domain": tune.lograndint(lower=2, upper=2**7 + 1),
                "init_value": 20,
            },
            "learning_rate": {
                "domain": tune.loguniform(lower=1 / 1024, upper=1.0),
                "init_value": 0.1,
            },
            "log_max_bin": {  # log transformed with base 2
                "domain": tune.lograndint(lower=3, upper=11),
                "init_value": 8,
            },
            "colsample_bytree": {
                "domain": tune.uniform(lower=0.01, upper=1.0),
                "init_value": 1.0,
            },
            "reg_alpha": {
                "domain": tune.loguniform(lower=1 / 1024, upper=1024),
                "init_value": 1 / 1024,
            },
            "reg_lambda": {
                "domain": tune.loguniform(lower=1 / 1024, upper=1024),
                "init_value": 1.0,
            },
        }


class LGBMEstimatorLimitDepth(LGBMEstimator):
    """LGBMEstimator with limited depth and higher floor for n_estimators."""

    @classmethod
    def search_space(cls, data_size, **params):
        return {
            "n_estimators": {
                "domain": tune.randint(lower=500, upper=5000),
                "init_value": 1000,
            },
            "num_leaves": {
                "domain": tune.randint(lower=6, upper=64),
                "init_value": 31,
            },
            "max_depth": {
                "domain": tune.randint(lower=4, upper=64),
                "init_value": 10,
            },
            "min_child_samples": {
                "domain": tune.lograndint(lower=2, upper=2**7 + 1),
                "init_value": 20,
            },
            "learning_rate": {
                "domain": tune.loguniform(lower=1 / 1024, upper=1.0),
                "init_value": 0.1,
            },
            "log_max_bin": {  # log transformed with base 2
                "domain": tune.lograndint(lower=3, upper=11),
                "init_value": 8,
            },
            "colsample_bytree": {
                "domain": tune.uniform(lower=0.01, upper=1.0),
                "init_value": 1.0,
            },
            "reg_alpha": {
                "domain": tune.loguniform(lower=1 / 1024, upper=1024),
                "init_value": 1 / 1024,
            },
            "reg_lambda": {
                "domain": tune.loguniform(lower=1 / 1024, upper=1024),
                "init_value": 1.0,
            },
        }


# class LassoSGD(SKLearnEstimator):
# 	"""Lasso with stochastic gradient descent."""

# 	def __init__(self, task="regression", n_jobs=None, **config):
# 		super().__init__(
# 			task,
# 			loss='squared_error',
# 			penalty='l1',
# 			**config
# 		)

# 		if task not in ["regression"]:
# 			raise ValueError("LassoSGD model is only for regression tasks")

# 		self.estimator_class = linear_model.SGDRegressor

# 	@classmethod
# 	def search_space(cls, data_size, task):
# 		space = {
# 			"alpha": {
# 				"domain": tune.loguniform(lower=1e-8, upper=1.0),
# 				"init_value": 1e-4,
# 			},
# 			"max_iter": {
# 				"domain": tune.lograndint(lower=500, upper=15000),
# 				"low_cost_init_value": 1000,
# 			},
# 			"tol": {
# 				"domain": tune.loguniform(lower=1e-7, upper=5e-3),
# 				"init_value": 1e-4,
# 			},
# 			"learning_rate": {
# 				"domain": tune.choice([
# 					'constant', 'invscaling', 'optimal', 'adaptive'
# 				]),
# 				"init_value": 'invscaling',
# 			},
# 			"eta0": {
# 				"domain": tune.loguniform(lower=1e-6, upper=0.1),
# 				"init_value": 0.01,
# 			},
# 			"power_t": {
# 				"domain": tune.uniform(lower=0.1, upper=0.5),
# 				"init_value": 0.25,
# 			},
# 			"n_iter_no_change": {
# 				"domain": tune.randint(3, upper=15),
# 				"init_value": 5,
# 			},
# 		}
# 		return space


# class LassoSGDEarlyStop(SKLearnEstimator):
# 	"""Lasso with stochastic gradient descent and early stopping."""

# 	def __init__(self, task="regression", n_jobs=None, **config):
# 		super().__init__(
# 			task,
# 			loss='squared_error',
# 			penalty='l1',
# 			early_stopping=True,
# 			**config
# 		)

# 		if task not in ["regression"]:
# 			raise ValueError("LassoSGD model is only for regression tasks")

# 		self.estimator_class = self._create_pipeline

# 	@staticmethod
# 	def _create_pipeline(scale=False, **kwargs):
# 		sgd_reg = linear_model.SGDRegressor(**kwargs)

# 		if scale is True:
# 			return pipeline.Pipeline([
# 				('scaler', StandardScaler()),
# 				('sgd', sgd_reg)
# 			])
# 		else:
# 			return sgd_reg

# 	@classmethod
# 	def search_space(cls, data_size, task):
# 		space = {
# 			"alpha": {
# 				"domain": tune.loguniform(lower=1e-12, upper=1e-6),
# 				"init_value": 1e-9,
# 			},
# 			"max_iter": {
# 				"domain": tune.lograndint(lower=500, upper=15000),
# 				"low_cost_init_value": 1000,
# 			},
# 			"tol": {
# 				"domain": tune.loguniform(lower=1e-7, upper=5e-3),
# 				"init_value": 1e-4,
# 			},
# 			"learning_rate": {
# 				"domain": tune.choice([
# 					'constant', 'invscaling', 'optimal', 'adaptive'
# 				]),
# 				"init_value": 'invscaling',
# 			},
# 			"eta0": {
# 				"domain": tune.loguniform(lower=1e-6, upper=0.1),
# 				"init_value": 0.01,
# 			},
# 			"power_t": {
# 				"domain": tune.uniform(lower=0.1, upper=0.5),
# 				"init_value": 0.25,
# 			},
# 			"n_iter_no_change": {
# 				"domain": tune.randint(3, upper=15),
# 				"init_value": 5,
# 			},
# 			"validation_fraction": {
# 				"domain": tune.uniform(lower=0.05, upper=0.25),
# 				"init_value": 0.1,
# 			},
# 			"scale": {
# 				"domain": tune.choice([
# 					True, False
# 				]),
# 				"init_value": True,
# 			},
# 		}
# 		return space


# class LassoSGDHuber(SKLearnEstimator):
# 	"""Lasso with stochastic gradient descent and Huber loss."""

# 	def __init__(self, task="regression", n_jobs=None, **config):
# 		super().__init__(
# 			task,
# 			loss='huber',
# 			penalty='l1',
# 			**config
# 		)

# 		if task not in ["regression"]:
# 			raise ValueError("LassoSGD model is only for regression tasks")

# 		self.estimator_class = linear_model.SGDRegressor

# 	@classmethod
# 	def search_space(cls, data_size, task):
# 		space = {
# 			"alpha": {
# 				"domain": tune.loguniform(lower=1e-8, upper=1.0),
# 				"init_value": 1e-4,
# 			},
# 			"max_iter": {
# 				"domain": tune.lograndint(lower=500, upper=15000),
# 				"low_cost_init_value": 1000,
# 			},
# 			"tol": {
# 				"domain": tune.loguniform(lower=1e-7, upper=5e-3),
# 				"init_value": 1e-4,
# 			},
# 			"learning_rate": {
# 				"domain": tune.choice([
# 					'constant', 'invscaling', 'optimal', 'adaptive'
# 				]),
# 				"init_value": 'invscaling',
# 			},
# 			"eta0": {
# 				"domain": tune.loguniform(lower=1e-6, upper=0.1),
# 				"init_value": 0.01,
# 			},
# 			"power_t": {
# 				"domain": tune.uniform(lower=0.1, upper=0.5),
# 				"init_value": 0.25,
# 			},
# 			"n_iter_no_change": {
# 				"domain": tune.randint(3, upper=15),
# 				"init_value": 5,
# 			},
# 			"epsilon": {
# 				"domain": tune.loguniform(lower=0.0001, upper=0.25),
# 				"init_value": 0.05,
# 			},
# 		}
# 		return space


# class LassoSGDHuberEarlyStop(SKLearnEstimator):
# 	"""Lasso with stochastic gradient descent, Huber loss, and early stopping."""

# 	def __init__(self, task="regression", n_jobs=None, **config):
# 		super().__init__(
# 			task,
# 			loss='huber',
# 			penalty='l1',
# 			early_stopping=True,
# 			**config
# 		)

# 		if task not in ["regression"]:
# 			raise ValueError("LassoSGD model is only for regression tasks")

# 		self.estimator_class = self._create_pipeline

# 	@staticmethod
# 	def _create_pipeline(scale=False, **kwargs):
# 		sgd_reg = linear_model.SGDRegressor(**kwargs)

# 		if scale is True:
# 			return pipeline.Pipeline([
# 				('scaler', StandardScaler()),
# 				('sgd', sgd_reg)
# 			])
# 		else:
# 			return sgd_reg

# 	@classmethod
# 	def search_space(cls, data_size, task):
# 		space = {
# 			"alpha": {
# 				"domain": tune.loguniform(lower=1e-12, upper=1e-6),
# 				"init_value": 1e-9,
# 			},
# 			"max_iter": {
# 				"domain": tune.lograndint(lower=500, upper=15000),
# 				"low_cost_init_value": 1000,
# 			},
# 			"tol": {
# 				"domain": tune.loguniform(lower=1e-7, upper=5e-3),
# 				"init_value": 1e-4,
# 			},
# 			"learning_rate": {
# 				"domain": tune.choice([
# 					'constant', 'invscaling', 'optimal', 'adaptive'
# 				]),
# 				"init_value": 'invscaling',
# 			},
# 			"eta0": {
# 				"domain": tune.loguniform(lower=1e-6, upper=0.1),
# 				"init_value": 0.01,
# 			},
# 			"power_t": {
# 				"domain": tune.uniform(lower=0.1, upper=0.5),
# 				"init_value": 0.25,
# 			},
# 			"n_iter_no_change": {
# 				"domain": tune.randint(3, upper=15),
# 				"init_value": 5,
# 			},
# 			"epsilon": {
# 				"domain": tune.loguniform(lower=0.0001, upper=0.25),
# 				"init_value": 0.05,
# 			},
# 			"validation_fraction": {
# 				"domain": tune.uniform(lower=0.05, upper=0.25),
# 				"init_value": 0.1,
# 			},
# 			"scale": {
# 				"domain": tune.choice([
# 					True, False
# 				]),
# 				"init_value": True,
# 			},
# 		}
# 		return space
