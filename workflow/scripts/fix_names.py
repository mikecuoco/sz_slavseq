#!/usr/bin/env python
__author__ = 'Michael Cuoco'

import pandas as pd
import sys, gc, traceback

def main():
	# read in the chromosome map
	chrom_map = pd.read_csv("resources/hs37d5_map.tsv", sep="\t", names=["hs37d5", "ann"])
	insert = pd.read_csv(snakemake.input[0], sep="\t")

	for name in chrom_map["ann"].to_list():
		insert.loc[insert["chr"] == name, "chr"] = chrom_map.loc[ chrom_map["ann"] == name, "hs37d5"].values[0]

	insert.to_csv(snakemake.output[0])

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

