#!/usr/bin/env python
# Created on: Aug 27, 2023 at 12:45:08 PM
__author__ = "Michael Cuoco"

import pandas as pd
import pyranges as pr


def remove_unused_categories(data):
    """Remove unused categories of columns Chromosome, donor_id, tissue_id, and cell_id"""

    for c in ["Chromosome", "donor_id", "tissue_id", "cell_id"]:
        data[c] = data[c].cat.remove_unused_categories()
    return data


def count_loci(data: pd.DataFrame, donor_knrgls: dict[str, pr.PyRanges]) -> int:
    "Count the number of loci in a dataset"

    # make cell_id a category, speed up grouping
    if data["cell_id"].dtype.name != "category":
        data.loc[:, "cell_id"] = data["cell_id"].astype("category")
    else:
        data.loc[:, "cell_id"] = data["cell_id"].cat.remove_unused_categories()

    # make donor_id a string if its not a string
    if data["donor_id"].dtype.name != "string":
        data.loc[:, "donor_id"] = data["donor_id"].astype(str)

    loci = 0
    for _, df in data.groupby(["cell_id"]):
        donor_id = df.donor_id.unique()[0]
        d = pr.PyRanges(df, int64=True)
        loci += len(donor_knrgls[donor_id].overlap(d).df)

    return loci
