import pandas as pd
from Bio.Seq import Seq
from pathlib import Path
from snakemake.utils import validate


validate(config, schema="../schemas/config.schema.yaml")

# read sample sheet
samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample": str, "donor": str})
    .set_index("sample", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")

# setup input/output for folds rule
def get_folds_input_samples(wildcards):
    my_samples = samples.loc[
        (samples["dna_type"] == wildcards.dna_type)
        & (samples["donor"] == wildcards.donor)
    ]["sample"]
    return expand(
        "results/flank_features/{donor}/{dna_type}/{sample}.pickle.gz",
        donor=wildcards.donor,
        dna_type=wildcards.dna_type,
        sample=my_samples,
    )


# get file of non-reference germline L1s
# if not from eul1db, should be a csv file with 3 or 4 columns:
# chrom, start, end, in_NRdb (optional)
def get_non_ref_l1(wildcards):
    db = config["ref"]["database"]
    NR_l1 = f"resources/{db}/windows.csv"

    if db != "eul1db":
        NR_df = pd.read_csv(NR_l1[0])
        validate(NR_df, schema="../schemas/non_ref_l1.schema.yaml")

    return NR_l1

num_folds = config["model"]["num_folds"]
fold_dirs = [f"fold_{fold}" for fold in range(num_folds)]
