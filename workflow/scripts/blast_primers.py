#!/usr/bin/env python
# Created on: Oct 28, 2023 at 12:49:30 AM
__author__ = "Michael Cuoco"


import sys, logging

# redirect stderr to log file
sys.stderr = open(snakemake.log[0], "w")  # type: ignore
logging.basicConfig(stream=sys.stderr, level=logging.INFO)

from subprocess import Popen
from pathlib import Path

# index the genome
db = Path(snakemake.output[0]).parent / snakemake.wildcards.genome

cmd = f"makeblastdb -in {snakemake.input.ref_fa} -dbtype 'nucl' -blastdb_version 5 -parse_seqids -out {db} >> {snakemake.log[0]} 2>&1"
logging.info(f"Running makeblastdb with command: {cmd}")
Popen(cmd, shell=True).wait()

# blast the primers
from tempfile import NamedTemporaryFile
import pandas as pd
import pyranges as pr

outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcov"

# add log statement
with NamedTemporaryFile() as tmp:
    cmd = f"blastn -query {snakemake.input.primer_fa} -task blastn -ungapped -no_greedy -db {db} -outfmt '{outfmt}' -out {tmp.name} >> {snakemake.log[0]} 2>&1"
    logging.info(f"Running blastn with command: {cmd}")
    Popen(cmd, shell=True).communicate()
    # if Popen(cmd, shell=True).wait() != 0:
    # 	raise SystemExit(f"Error: blastn exited with non-zero exit code: {exit}")
    df = pd.read_csv(tmp.name, sep="\t", header=None, names=outfmt.split(" ")[1:])


def switch_start_end(x):
    if x.Start > x.End:
        x.Start, x.End = x.End, x.Start
    return x


new_names = {
    "qseqid": "Name",
    "sseqid": "Chromosome",
    "sstart": "Start",
    "send": "End",
    "sstrand": "Strand",
    "bitscore": "Score",
}

df = df[list(new_names)].rename(new_names, axis=1).apply(switch_start_end, axis=1)
df["Strand"] = df["Strand"].str.replace("minus", "-")
df["Strand"] = df["Strand"].str.replace("plus", "+")

df_capture = pr.PyRanges(df[df["Name"] == "L1_capture_probe"])
df_pcr = pr.PyRanges(df[df["Name"] == "L1_PCR_primer"])

df = df_capture.extend(1000).overlap(df_pcr).df
df["Name"] = "Predicted_SLAVseq_fragment"

# save as bed
logging.info(f"Saving blast results to {snakemake.output.bed}")
pr.PyRanges(df).to_bed(snakemake.output.bed)  # type: ignore

sys.stderr.close()
