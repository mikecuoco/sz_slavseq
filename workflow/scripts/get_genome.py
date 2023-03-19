#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import sys, tempfile
from pathlib import Path
from snakemake.shell import shell


sys.stderr = open(snakemake.log[0], "w")

# check if file is gzipped, if so, unzip to temp file
tmp = tempfile.NamedTemporaryFile(suffix=".fa")
if Path(snakemake.input[0]).suffix == ".gz":
    shell("gunzip -c {snakemake.input[0]} > {tmp.name}")
    snakemake.input[0] = tmp.name

# get region(s) if specified
if snakemake.config["genome"]["region"] != "all":

    if isinstance(snakemake.config["genome"]["region"], list):
        region = " ".join(snakemake.config["genome"]["region"])
    else:
        region = snakemake.config["genome"]["region"]

    print("Extracting region(s): {} from {}".format(region, snakemake.input[0]))
    shell("samtools faidx {snakemake.input[0]} {region} > {snakemake.output.fa}")

else:
    shell("cp {snakemake.input[0]} {snakemake.output.fa}")

# index
shell("samtools faidx {snakemake.output.fa}")

tmp.close()
sys.stderr.close()
