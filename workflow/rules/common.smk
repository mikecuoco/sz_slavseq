import pandas as pd 
from Bio.Seq import Seq

# read sample sheet 
samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

# TODO add sample sheet validation schema
# validate(samples, schema="../schemas/samples.schema.yaml")

num_folds=config["model"]["num_folds"]
b_names = ["X_train.pickle.gz", "X_test.pickle.gz", "Y_train.pickle.gz", "Y_test.pickle.gz"]
fold_files = []
for fold in range(num_folds):
    for f in b_names:
        path = f"fold_{fold}/" + f
        fold_files.append(path)
