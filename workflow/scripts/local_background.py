#!/usr/bin/env python
# Created on: Nov 19, 2023 at 4:21:51 PM
__author__ = "Michael Cuoco"

import logging

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import pyarrow.parquet as pq
import pandas as pd
from pysam import AlignmentFile
from pyslavseq.sliding_window import SlidingWindow, read_to_namedtuple
from pyslavseq.features import features

logger.info(f"Generating local background features from {snakemake.input.bam}")  # type: ignore
regions = pq.read_table(snakemake.input["pqt"]).to_pandas()  # type: ignore

bg = [5e3, 1e4, 2e4]
side = ["left", "right"]
bg_dfs = {(bg, s): pd.DataFrame() for b in bg for s in side}

with AlignmentFile(snakemake.input["bam"], "rb") as bam:  # type: ignore
    sw = SlidingWindow(bam)
    for (s, f), df in bg_dfs.items():
        for r in regions.itertuples():
            # setup window
            w = {"Chromosome": r.Chromosome}
            w["Start"] = r.Start - 5000 if flank == "left" else r.End
            w["End"] = r.Start if flank == "left" else r.End + 5000

            # get reads
            reads = bam.fetch(
                w["Chromosome"], w["Start"], w["End"]
            )  # this is faster the a sliding window
            reads = filter(sw.read_filter, reads)
            reads = filter(lambda x: x.is_read1, reads)
            w["reads"] = list(map(read_to_namedtuple, reads))

            # get features, add to df
            f = features(w)
            df.loc[r.Index, f.keys()] = f

for flank, df in bg_dfs.items():
    df.drop(columns=["Chromosome", "Start", "End"], inplace=True)
    df.columns = [f"{flank}_{c}" for c in df.columns]
    regions = regions.join(df)

regions.to_parquet(snakemake.output[0])  # type: ignore
