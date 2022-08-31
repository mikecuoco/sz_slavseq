#!/usr/bin/env python
# adapted from Ricardo's jupyter notebook from rf pipeline: pipeline/eul1db/to_windows/_h/main.ipynb
__author__ = 'Michael Cuoco'

# TODO: add comments 

import pandas as pd
from pyslavseq.genome import Interval, Genome
import pyslavseq.genome.interval_generator as ig
import snakemake as sm
import logging

def download_rmsk(ref):
    if "38" in ref:
        url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_rm.out.gz"
        sm.shell(f"wget -O- --no-config -q {url} > resources/{snakemake.wildcards.ref}rmsk.out.gz")
        sm.shell("perl -e 'use Data::Dumper; print Dumper(\@INC)'") # print value of @INC
        sm.shell(f"rmToUCSCTables.pl -out resources/{snakemake.wildcards.ref}/rmsk.out.gz") # utility from RepeatMasker
        sm.shell(f"gzip -c resources/{snakemake.wildcards.ref}/rmsk.out.tsv > {snakemake.output.rmsk}")
    elif "37" in ref:
        url="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz"
        sm.shell(f"wget -O- --no-config -q {url} > {snakemake.output.rmsk}")


def read_rmsk():

    # read the rmsk file
    df0 = pd.read_csv(snakemake.output.rmsk, compression="gzip", sep="\t", header = None, usecols = [5,6,7,9,10])

    # filter for rep_names
    rep_names = ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]
    logging.info(f"Filtering for rep_names: {rep_names}")
    df0 =  df0.loc[df0.iloc[:,4].isin(rep_names)]

    # save to new dataframe
    df1 = pd.DataFrame()
    df1['chrom'] = df0.iloc[:,0]
    # set start positions depending on strand
    df1['pos'] = df0.apply(lambda x: x.iloc[2] if x.iloc[3] == '+' else x.iloc[1], axis=1)
    df1['l1_name'] = df0.iloc[:,4]

    return df1


def main():

    # setup logging
    logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

    download_rmsk(snakemake.wildcards.ref)
    df = read_rmsk()
    l1_pos = set()

    for (_, chrom, pos) in df[['chrom', 'pos']].itertuples():
        l1_pos.update([Interval(chrom, pos, pos)])
        
    genome = Genome(snakemake.input[0])
    xx = list(ig.windows_overlapping_intervals(genome, l1_pos, 750, 250))

    refl1 = pd.DataFrame.from_records((x.as_tuple() for x in xx), columns=['chrom', 'start', 'end']).set_index(['chrom', 'start', 'end'])
    refl1['reference_l1hs_l1pa2_6'] = True

    refl1.to_csv(snakemake.output.ref_l1)

if __name__ == '__main__':
    main()