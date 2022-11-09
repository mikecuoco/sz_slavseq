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


# handle specified region
region = (
    "".join(config["genome"]["region"])
    if isinstance(config["genome"]["region"], list)
    else config["genome"]["region"]
)
region_name = f"_{region}" if region != "all" else ""

# make name of bed file for get_eul1db rule
def get_liftover_input(wildcards):
    if wildcards.db == "eul1db":
        return "resources/hg19/hg19_eul1db_insertions.bed"
    elif wildcards.db == "dbVar":
        return "resources/hs38DH/hs38DH_dbVar_insertions.bed"


def get_fixnames_input(wildcards):
    if wildcards.ref == "hs37d5":
        return f"resources/hg19/hg19_{wildcards.db}_insertions.bed"
    else:
        return (
            f"resources/{wildcards.ref}/{wildcards.ref}_{wildcards.db}_insertions.bed"
        )
