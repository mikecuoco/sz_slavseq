#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys
import pandas as pd
import pickle
from sklearn.utils import shuffle
from sklearn.metrics import precision_recall_curve, auc, confusion_matrix


sys.stderr = open(snakemake.log[0], "w")

# get label encoder
with open(snakemake.input.label_encoder, "rb") as f:
    le = pickle.load(f)

# get predictions
with open(snakemake.input.proba, "rb") as f:
    proba = pickle.load(f)

with open(snakemake.input.pred, "rb") as f:
    pred = pickle.load(f)

# get labels
with open(snakemake.input.labels, "rb") as f:
    y = pickle.load(f)

# get PR curves and confusion matrices in a loop
pr_list = []; cm_dict = {}
for stage in ["train", "test", "test_shuffled"]:
    cm_dict[stage] = {}
    for fold in proba.keys():
        if stage == "test_shuffled":
            cm_dict[stage][fold] = confusion_matrix(
                y_true = shuffle(le.inverse_transform(y[fold]["test"])),
                y_pred = le.inverse_transform(pred[fold]["test"]), 
                labels = le.classes_
            )
        else:
            cm_dict[stage][fold] = confusion_matrix(
                y_true = le.inverse_transform(y[fold][stage]), 
                y_pred = le.inverse_transform(pred[fold][stage]), 
                labels = le.classes_
            )
        for label in ["KNRGL", "RL1", "OTHER"]:

            pr = pd.DataFrame()

            if stage == "test_shuffled":
                pr["precision"], pr["recall"], _ = precision_recall_curve(
                    shuffle(y[fold]["test"]),
                    proba[fold]["test"][:, le.transform([label])[0]].round(2),
                    pos_label=le.transform([label])[0],
                )
            else:
                pr["precision"], pr["recall"], _ = precision_recall_curve(
                    y[fold][stage],
                    proba[fold][stage][:, le.transform([label])[0]].round(2),
                    pos_label=le.transform([label])[0],
                )
                
            pr["stage"], pr["label"], pr["fold"] = stage, label, fold + 1
            pr["model_id"] = snakemake.wildcards.model_id
            pr["auprc"] = auc(pr["recall"], pr["precision"])
            pr_list.append(pr)

# save confusion matrices
with open(snakemake.output.confusion, "wb") as f:
    pickle.dump(cm_dict, f)

# save PR curves
pr = pd.concat(pr_list).reset_index(drop=True)
pr.to_pickle(snakemake.output.prcurve)


sys.stderr.close()
