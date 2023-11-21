#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from numpy import max
import pyBigWig

# get first cell to set schema
cell = pq.read_table(snakemake.input.cells[0]).to_pandas()
cell["cell_id"] = Path(snakemake.input.cells[0]).name.rstrip(".pqt")
schema = pa.Schema.from_pandas(cell)

# add en_motif scores to schema
schema = schema.append(pa.field("en_pos_score", pa.float32()))
schema = schema.append(pa.field("en_neg_score", pa.float32()))

# open bigwigs
en_pos_bw = pyBigWig.open(snakemake.input.en_motif_pos, "r")
en_neg_bw = pyBigWig.open(snakemake.input.en_motif_neg, "r")

# stream this to reduce memory usage
logger.info(f"Reading regions for donor {snakemake.wildcards.donor}..")
try:
    with pq.ParquetWriter(snakemake.output[0], schema, compression="gzip") as writer:
        for f in snakemake.input.cells:
            logger.info(f"Reading {f}")
            cell = pq.read_table(f).to_pandas().query("n_ref_reads == 0")
            cell = cell.loc[cell["Chromosome"].isin([f"chr{i}" for i in range(1, 23)])]
            cell["cell_id"] = Path(f).name.rstrip(".pqt")

            cell["en_pos_score"] = cell.apply(
                lambda x: max(en_pos_bw.values(x["Chromosome"], x["Start"], x["End"])),
                axis=1,
            )
            cell["en_neg_score"] = cell.apply(
                lambda x: max(en_neg_bw.values(x["Chromosome"], x["Start"], x["End"])),
                axis=1,
            )

            writer.write_table(pa.Table.from_pandas(cell, schema=schema))
except:
    logger.error("Error/Interrupt detected, deleting output file")
    Path(snakemake.output[0]).unlink()
    raise Exception("Error/Interrupt detected, deleting output file")
