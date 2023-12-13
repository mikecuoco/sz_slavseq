#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import logging

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from numpy import max
import pyBigWig
import pyranges as pr
from pyslavseq.preprocessing import get_stratification

# get first cell to set schema
cell = pq.read_table(snakemake.input.cells[0]).to_pandas()  # type: ignore
cell["cell_id"] = Path(snakemake.input.cells[0]).name.rstrip(".pqt")  # type: ignore
schema = pa.Schema.from_pandas(cell)

# add en_motif scores to schema
schema = schema.append(pa.field("en_pos_score", pa.float32()))
schema = schema.append(pa.field("en_neg_score", pa.float32()))

# open bigwigs
en_pos_bw = pyBigWig.open(snakemake.input.en_motif_pos, "r")  # type: ignore
en_neg_bw = pyBigWig.open(snakemake.input.en_motif_neg, "r")  # type: ignore

logger.info(f"Reading non-unique mappable regions for {snakemake.wildcards.genome}...")  # type: ignore
if "chm" in snakemake.wildcards.genome:  # type: ignore
    nonunique = get_stratification("chm13", "nonunique")
else:
    nonunique = get_stratification("grch38", "nonunique")
nonunique = pr.PyRanges(nonunique)

# stream this to reduce memory usage
logger.info(f"Reading regions for donor {snakemake.wildcards.donor}..")  # type: ignore
try:
    with pq.ParquetWriter(snakemake.output[0], schema, compression="gzip") as writer:  # type: ignore
        for f in snakemake.input.cells:  # type: ignore
            logger.info(f"Reading {f}")
            cell = (
                pq.read_table(f)
                .to_pandas()
                .query("max_mapq > 0")
                .reset_index(drop=True)
            )

            # keep autosomes
            cell = cell.loc[cell["Chromosome"].isin([f"chr{i}" for i in range(1, 23)])]

            # remove non-unique mappable regions
            before = cell.shape[0]
            cell = (
                pr.PyRanges(cell).overlap(nonunique, how="containment", invert=True).df
            )
            logger.info(
                f"Removed {before-cell.shape[0]} regions contained in non-unique mappable regions"
            )

            # add cell id
            cell["cell_id"] = Path(f).name.rstrip(".pqt")

            # add en_motif scores
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
    Path(snakemake.output[0]).unlink()  # type: ignore
    raise Exception("Error/Interrupt detected, deleting output file")
