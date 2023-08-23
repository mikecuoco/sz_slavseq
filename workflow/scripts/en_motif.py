#!/usr/bin/env python
# Created on: Aug 1, 2023 at 12:12:15 PM
# followed BioPython motif guide: https://biopython-tutorial.readthedocs.io/en/latest/notebooks/14%20-%20Sequence%20motif%20analysis%20using%20Bio.motifs.html#
__author__ = "Michael Cuoco"

from Bio import motifs, SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from collections import defaultdict
import sys
import pyarrow as pa
import pyarrow.parquet as pq

sys.stderr = open(snakemake.log[0], "w")

## CREATE MOTIF
# Engineered L1 Insertion Coordinates, Characteristics, and Sequences from Flasch et al. 2019 PMID: 30955886
FLASCH_2019_SUPP1 = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6558663/bin/NIHMS1523125-supplement-8.xls"
flasch = pd.read_excel(FLASCH_2019_SUPP1)
instances = [Seq(site.replace("/", "")) for site in flasch["L1 EN Cleavage"].values]
m = motifs.create(instances)

# compute genomic background frequencies
bg = defaultdict(list)
with open(snakemake.input[0], "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        for base in "ACGT":
            bg[base].append(record.seq.count(base) / len(record.seq))

# average across all chromosomes
for base in bg.keys():
    bg[base] = sum(bg[base]) / len(bg[base])

# create PWM and PSSM
pwm = m.counts.normalize()
pssm = pwm.log_odds(bg)
rpssm = pssm.reverse_complement()

## SCORE GENOME
schema = pa.schema(
    [
        pa.field("Chromosome", pa.string()),
        pa.field("Start", pa.int64()),
        pa.field("End", pa.int64()),
        pa.field("pos_score_q0", pa.float32()),
        pa.field("pos_score_q0.25", pa.float32()),
        pa.field("pos_score_q0.5", pa.float32()),
        pa.field("pos_score_q0.75", pa.float32()),
        pa.field("pos_score_q1", pa.float32()),
        pa.field("pos_score_mean", pa.float32()),
        pa.field("neg_score_q0", pa.float32()),
        pa.field("neg_score_q0.25", pa.float32()),
        pa.field("neg_score_q0.5", pa.float32()),
        pa.field("neg_score_q0.75", pa.float32()),
        pa.field("neg_score_q1", pa.float32()),
        pa.field("neg_score_mean", pa.float32()),
    ]
)

with open(snakemake.input[0], "r") as infile, pq.ParquetWriter(
    snakemake.output[0], schema, compression="gzip"
) as writer:
    for record in SeqIO.parse(infile, "fasta"):
        if record.id not in [f"chr{c}" for c in range(1, 23)]:
            continue

        # calculate scores for each base this chromosome
        scores = pd.DataFrame(
            {
                "Chromosome": record.id,
                "pos_score": pssm.calculate(record.seq),
                "neg_score": rpssm.calculate(record.seq),
            }
        )
        scores["Start"] = scores.index + 1
        scores["End"] = scores["Start"] + pssm.length
        min_score = min(scores["pos_score"].min(), scores["neg_score"].min())
        scores.fillna(min_score - 100, inplace=True)

        # generate windows for this chromosome
        reflen = len(record.seq)
        step = 250
        size = 750
        window_scores = []
        for start in range(0, reflen + 1, step):
            end = start + size if start + size < reflen else reflen
            out = {
                "Chromosome": record.id,
                "Start": start,
                "End": end,
            }

            for strand in ["pos_score", "neg_score"]:
                try:
                    quantiles = np.quantile(
                        scores.loc[start : (end - len(pssm)), strand],
                        [0, 0.25, 0.5, 0.75, 1],
                    )
                except:
                    IndexError
                    import pdb

                    pdb.set_trace()

                for n, q in zip([0, 0.25, 0.5, 0.75, 1], quantiles):
                    out[strand + "_q" + str(n)] = q
                out[strand + "_mean"] = scores.loc[
                    start : (end - len(pssm)), strand
                ].mean()

            window_scores.append(out)

        writer.write_table(
            pa.Table.from_pandas(pd.DataFrame(window_scores), schema=schema)
        )

        # write to file

sys.stderr.close()
