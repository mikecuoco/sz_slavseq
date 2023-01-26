from time import time
import sys
import numpy as np
import pandas as pd

class Model:
	def __init__(self, df: pd.DataFrame, features: list, min_reads: int):
		self.df_ = df[
			df["all_reads.count"] >= min_reads
		].reset_index()  # remove windows with too few reads
		self.min_reads_ = min_reads
		self.features_ = features
		self.classes_ = self.df_["label"].unique()

	def split(
		self,
		kf,
		split_by: str,
		test_sampler=None,
	):

		# set configuration
		self.kf_ = kf
		self.test_sampler_ = test_sampler
		self.split_by_ = split_by

		# split the data into folds
		groups = self.df_[split_by]

		self.folds_ = {}
		for fold, (train_index, test_index) in enumerate(
			self.kf_.split(self.df_, self.df_["label"], groups=groups)
		):
			print(f"Generating fold {fold+1}...")
			self.folds_[fold] = {}

			# make train set
			self.folds_[fold]["train"] = self.df_.loc[train_index, :]
			assert self.folds_[fold]["train"]["label"].nunique() == len(
				self.classes_
			), "All labels are not present in the train set."
			print("Train set size:", file=sys.stderr)
			print(self.folds_[fold]['train']["label"].value_counts(), file=sys.stderr)

			# make test set
			if test_sampler:
				self.folds_[fold]["test"], _ = test_sampler.fit_resample(
					self.df_.loc[test_index, :], self.df_.loc[test_index, "label"]
				)
			else:
				self.folds_[fold]["test"] = self.df_.loc[test_index, :]
			assert self.folds_[fold]["test"]["label"].nunique() == len(
				self.classes_
			), "All labels are not present in the test set."
			print("Test set size: ", file=sys.stderr)
			print(self.folds_[fold]['test']["label"].value_counts(), file=sys.stderr)

		delattr(self, "df_") # save space

		return self

	def train_test(self, clf: object, train_sampler=None, id: str = None):

		# set configuration
		self.train_sampler_ = train_sampler
		self.id_ = id

		for fold in self.folds_:

			# resample the training set if specified
			if train_sampler:
				self.folds_[fold]["train"], _ = train_sampler.fit_resample(
					self.folds_[fold]["train"], self.folds_[fold]["train"]["label"]
				)
				print("Resampled train set size: ", file=sys.stderr)
				print(self.folds_[fold]['test']["label"].value_counts(), file=sys.stderr)
				assert self.folds_[fold]["train"]["label"].nunique() == len(
					self.classes_
				), "All labels are not present in the training set."

			# fit the model
			print(f"Training on fold {fold+1}...", file=sys.stderr)
			st = time()
			clf.fit(
				self.folds_[fold]["train"][self.features_],
				self.folds_[fold]["train"]["label"],
			)
			sp = time()
			print(f"Finished in {sp-st:.2f} seconds\n")

			# save the model
			self.folds_[fold]["clf"] = clf

			for stage in ["train", "test"]:

				self.folds_[fold][stage]["pred"] = clf.predict(
					self.folds_[fold][stage].loc[:, self.features_]
				)

				proba = clf.predict_proba(self.folds_[fold][stage].loc[:, self.features_])

				for label in self.classes_:
					i, = np.where(clf.classes_ == label) # get index of label
					self.folds_[fold][stage][f"proba_{label}"] = proba[:, i]
			
		return self
