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
