import pandas as pd
from pathlib import Path
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

# read sample sheet
samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample": str, "donor": str})
    .set_index(["sample", "donor", "dna_type"], drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")
