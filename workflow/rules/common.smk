import pandas as pd
from Bio.Seq import Seq
from pathlib import Path
from snakemake.utils import validate

outdir = config["outdir"]
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
