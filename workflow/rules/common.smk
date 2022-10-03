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
gen_ref_basename = ref

# handle trimming for a region
region = config["genome"]["region"]
if region != "all":
    gen_ref_basename = f"{ref}_{region}"

# handle non-reference L1 conversion to bed
# TODO: move rule input to variable here, make script amenable to any input
db = config["germline_line1"]["source"]

# final desired bed file, including wildcards
non_ref_l1_bed_final = "resources/{ref}/{ref}_{db}_insertions.bed"

# make name of bed file for get_eul1db rule
if db == "eul1db":
    if ref == "hs37d5":
        non_ref_l1_bed = f"resources/{ref}/{ref}_{db}_insertions.bed"
    else:
        non_ref_l1_bed = f"resources/hg19/hg19_{db}_insertions.bed"

# get raw non-reference germline L1s file
# if not from eul1db, should be a csv file with 3 or 4 columns:
# chrom, start, end, in_NRdb (optional)
non_ref_l1_windows = f"resources/{ref}/{db}_windows.csv"
if db != "eul1db":
    NR_df = pd.read_csv(non_ref_l1_windows)
    validate(NR_df, schema="../schemas/non_ref_l1.schema.yaml")


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


# handle number of folds
num_folds = config["model"]["num_folds"]
fold_dirs = [f"fold_{fold}" for fold in range(num_folds)]
