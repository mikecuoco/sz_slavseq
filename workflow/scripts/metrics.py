#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = 'Michael Cuoco'

import pandas as pd
import pickle
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt
import seaborn as sns

def get_prc(num_folds, le, models, label_files, proba_files, outimg):
	# get data and plot in loop
	df = pd.DataFrame()
	for fold in range(1,num_folds+1):
		for stage in ["train", "test"]: 

			with open(label_files[stage].format(fold=fold), "rb") as f:
				y = pickle.load(f)
			y = le.inverse_transform(y)
			y[y != "KNRGL"] = "OTHER" # make binary comparison

			for model in models:
				_df = pd.DataFrame()
				
				with open(proba_files[stage].format(fold=fold, model=model), "rb") as f:
					y_proba = pickle.load(f)[:,1]

				_df["precision"], _df["recall"], _ = precision_recall_curve(le.transform(y), y_proba, pos_label=le.transform(["KNRGL"])[0])

				_df["stage"] = stage
				_df["model"] = model
				_df["fold"] = fold
				df = pd.concat([df, _df]).reset_index(drop=True)

	sns.set_style("ticks")
	fig = sns.relplot(data=df, x="recall", y="precision", hue="model", col="stage", kind="line", markers=True)
	fig.set(ylim = [0,1], xlim = [0,1])
	plt.savefig(outimg, dpi=300)

if __name__ == "__main__":

	# setup inputs
	with open(snakemake.input.label_encoder, "rb") as f:
		le = pickle.load(f)

	label_files = {"train": snakemake.input.train_labels, "test": snakemake.input.test_labels}
	proba_files = {"train": snakemake.input.train_probabilities, "test": snakemake.input.test_probabilities}

	get_prc(snakemake.params.num_folds, le, label_files, proba_files, snakemake.output.prcurve)