import pandas as pd 

# read sample sheet
samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_id": str})
    .set_index("sample_id", drop=False)
    .sort_index()
)

# TODO add sample sheet validation schema
# validate(samples, schema="../schemas/samples.schema.yaml")
