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

# setup input/output for folds rule
def get_folds_input_samples(wildcards):
    my_samples = samples.loc[(samples['dna_type'] == wildcards.dna_type) & (samples['donor'] == int(wildcards.donor))]['sample']
    return expand("results/flank_features/{donor}/{dna_type}/{sample}.pickle.gz", donor=wildcards.donor, dna_type=wildcards.dna_type, sample=my_samples)

num_folds=config["model"]["num_folds"]
fold_dirs = [f"fold_{fold}" for fold in range(num_folds)]
