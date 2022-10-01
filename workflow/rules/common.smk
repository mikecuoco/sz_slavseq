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


def get_cutadapt_input(wildcards):
    sample = samples.loc[wildcards.sample]

    if "R1" in sample:
        return [sample["R1"], sample["R2"]]
    else:
        accession = sample["sra"]
        return expand(
            "results/fastq/{accession}_{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )


# handle conditional alterations to reference genome
ref = config["genome"]["build"]
gen_ref_basename = "genome"

# handle trimming for a region
region = config["genome"]["region"]
if region != "all":
    gen_ref_basename = f"genome_{region}"

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


def get_non_ref_l1_for_liftover(wildcards):
    db = config["germline_line1"]["source"]

    if db == "eul1db":
        return f"resources/{db}/insertions_hs37d5.bed"


# get file of non-reference germline L1s
# if not from eul1db, should be a csv file with 3 or 4 columns:
# chrom, start, end, in_NRdb (optional)
def get_non_ref_l1(wildcards):
    db = config["germline_line1"]["source"]
    NR_l1 = f"resources/{db}/windows.csv"

    if db != "eul1db":
        NR_df = pd.read_csv(NR_l1[0])
        validate(NR_df, schema="../schemas/non_ref_l1.schema.yaml")

    return NR_l1


num_folds = config["model"]["num_folds"]
fold_dirs = [f"fold_{fold}" for fold in range(num_folds)]
