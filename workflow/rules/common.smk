import pandas as pd 
from Bio.Seq import Seq
from pathlib import Path

# read sample sheet 
samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

# TODO add sample sheet validation schema
# validate(samples, schema="../schemas/samples.schema.yaml")

# setup output folders for folds rule
num_folds=config["model"]["num_folds"]
fold_dirs = [f"fold_{fold}" for fold in range(num_folds)]
