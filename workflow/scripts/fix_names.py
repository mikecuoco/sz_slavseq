#!/usr/bin/env python
__author__ = 'Michael Cuoco'

from Bio import SeqIO
import pandas as pd
import pysam
import os
import logging
import snakemake as sm

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

if "37" in snakemake.wildcards.ref:
	# read in the chromosome map
	# TODO: make this a parameter specified in snakemake
	chrom_map = pd.read_csv("resources/hs37d5_map.tsv", sep="\t", names=["hs37d5", "ann"])
	insert = pd.read_csv(snakemake.input[0], sep="\t", header=True)

	for name in chrom_map["ann"].to_list():
		insert.loc[ insert["chr"].str.contains(name), "chr"] = chrom_map.loc[ chrom_map["ann"] == name, "hs37d5"].values[0]
		
else:
	os.rename(snakemake.input[0], snakemake.output[0])
