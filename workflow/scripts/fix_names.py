#!/usr/bin/env python
__author__ = 'Michael Cuoco'

from Bio import SeqIO
import pandas as pd
import pysam
import os
import logging

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

if "37" in snakemake.wildcards.ref:
	# read in the chromosome map
	chrom_map = pd.read_csv("workflow/scripts/hs37d5_map.tsv", sep="\t", header=None)

	# save filenames to objects
	original_file = snakemake.input.fa[0]
	corrected_file = snakemake.output.fa

	with open(corrected_file, 'w') as corrected:
		records = SeqIO.parse(original_file, 'fasta')
		for record in records:

			# check if the record is in the chromosome map and change it
			if chrom_map.iloc[:,0].str.contains(record.id).any():
				logging.info(f"found {record.id} in map")
				record.id = chrom_map.iloc[:,1][chrom_map.iloc[:,0] == record.id].values[0]
				logging.info(f"change to {record.id}")
			else:
				logging.info(f"did not find {record.id} in map")

			# write the record to the output file
			SeqIO.write(record, corrected, 'fasta')
else:
	os.rename(snakemake.input.fa[0], snakemake.output.fa)

pysam.faidx(snakemake.output.fa)
fai = pd.read_csv(snakemake.output.fai, sep="\t", header=None)
fai.iloc[:,0:1].to_csv(snakemake.output.chromsizes, sep="\t", header=False, index=False)