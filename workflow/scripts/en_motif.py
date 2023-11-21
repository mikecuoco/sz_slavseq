#!/usr/bin/env python
# Created on: Aug 1, 2023 at 12:12:15 PM
# followed BioPython motif guide: https://biopython-tutorial.readthedocs.io/en/latest/notebooks/14%20-%20Sequence%20motif%20analysis%20using%20Bio.motifs.html#
__author__ = "Michael Cuoco"

import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

from Bio import motifs, SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from collections import defaultdict
import pyBigWig


CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

## CREATE MOTIF
# Engineered L1 Insertion Coordinates, Characteristics, and Sequences from Flasch et al. 2019 PMID: 30955886
logger.info("Loading L1 EN cleavage sites from Flasch et al. 2019")
FLASCH_2019_SUPP1 = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6558663/bin/NIHMS1523125-supplement-8.xls"
flasch = pd.read_excel(FLASCH_2019_SUPP1)
instances = [Seq(site.replace("/", "")) for site in flasch["L1 EN Cleavage"].values]
m = motifs.create(instances)

# compute genomic background frequencies
logger.info(f"Computing genomic background frequencies from {snakemake.input[0]}")
bg = defaultdict(list)
with open(snakemake.input[0], "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        if record.id not in CHROMOSOMES:
            continue
        for base in "ACGT":
            bg[base].append(record.seq.count(base) / len(record.seq))

# average across all chromosomes
for base in bg.keys():
    bg[base] = sum(bg[base]) / len(bg[base])

# create PWM and PSSM
logger.info("Creating PWM and PSSM")
pwm = m.counts.normalize()
pssm = pwm.log_odds(bg)
rpssm = pssm.reverse_complement()

# setup output files
pos_wig = pyBigWig.open(snakemake.output.pos, "w")
neg_wig = pyBigWig.open(snakemake.output.neg, "w")

header = []
for record in SeqIO.parse(snakemake.input[0], "fasta"):
    if record.id not in CHROMOSOMES:
        continue
    header.append((record.id, len(record.seq)))

pos_wig.addHeader(header)
neg_wig.addHeader(header)

## SCORE GENOME
for record in SeqIO.parse(snakemake.input[0], "fasta"):
    if record.id not in CHROMOSOMES:
        continue

    # calculate scores for each base this chromosome
    pos_score = pssm.calculate(record.seq)
    neg_score = rpssm.calculate(record.seq)

    # write scores to wiggle
    logger.info(f"Writing scores for {record.id}")
    pos_wig.addEntries(record.id, 0, values=pos_score, span=pssm.length, step=1)
    neg_wig.addEntries(record.id, 0, values=neg_score, span=rpssm.length, step=1)

# close files
pos_wig.close()
neg_wig.close()
