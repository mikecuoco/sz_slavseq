#!/usr/bin/env python
# Created on: Nov 17, 2023 at 6:45:57 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

from pathlib import Path
import numpy as np
import pandas as pd
import pyranges as pr
import pyBigWig
from pyslavseq.preprocessing import collate_labels

if __name__ == "__main__":
    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )

    logger = logging.getLogger(__name__)

    regions = pd.read_parquet(snakemake.input.regions).set_index(
        ["Chromosome", "Start", "End"]
    )

    for a in map(Path, snakemake.input.annotations):  # type: ignore
        name = a.name.split(".")[1]
        # check if file is empty
        if a.stat().st_size == 0:
            regions[name] = False
            continue
        anno = pr.read_bed(str(a)).df.set_index(["Chromosome", "Start", "End"])
        regions[name] = regions.index.isin(anno.index)
    regions.reset_index(inplace=True)

    regions.to_parquet(snakemake.output[0])  # type: ignore
