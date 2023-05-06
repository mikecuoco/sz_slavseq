import pandas as pd
from pathlib import Path
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

# read sample sheets
samples = pd.read_csv(
    config["samples"],
    sep="\t",
    dtype={"sample_id": str, "tissue_id": str, "donor_id": str},
)
validate(samples, schema="../schemas/samples.schema.yaml")
donors = pd.read_csv(config["donors"], sep="\t", dtype={"donor_id": str})
validate(donors, schema="../schemas/donors.schema.yaml")

# merge sample sheets
samples = samples.merge(donors, on=["donor_id"]).set_index(
    ["donor_id", "sample_id"], drop=False
)

# create donor sheet
donors = donors.set_index("donor_id", drop=False)
