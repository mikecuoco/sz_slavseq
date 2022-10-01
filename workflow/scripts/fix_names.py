#!/usr/bin/env python
__author__ = 'Michael Cuoco'

from Bio import SeqIO
import pandas as pd
import pysam
import sys, gc, traceback
import snakemake as sm
import logging 

def main():

	logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)
	# read in the chromosome map
	chrom_map = pd.read_csv("resources/hs37d5_map.tsv", sep="\t", header=None)

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


	# samtools faidx
	pysam.faidx(snakemake.output.fa)
	sm.shell("cut -f 1,2 {snakemake.output.fai} > {snakemake.output.chromsizes}")

if __name__ == '__main__':
    try:
        main()

    except:  # catch *all* exceptions
        sys.stderr = open(snakemake.log[0], 'w')
        traceback.print_exc()
        sys.stderr.close()

    finally:
        # cleanup code in here
        gc.collect()