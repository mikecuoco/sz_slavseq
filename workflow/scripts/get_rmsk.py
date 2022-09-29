#!/usr/bin/env python
# adapted from Ricardo's jupyter notebook from rf pipeline: pipeline/eul1db/to_windows/_h/main.ipynb
__author__ = 'Michael Cuoco, Rohini Gadde'

import pandas as pd
from pyslavseq.genome import Interval, Genome
import pyslavseq.genome.interval_generator as ig
import snakemake as sm
import logging
import sys, gc, traceback

def read_rmsk():
    # read the rmsk file
    df0 = pd.read_csv(snakemake.input['rmsk'], skiprows=3, delim_whitespace=True, 
        names=["chr", "start", "end", "strand", "repeat"], usecols=[4,5,6,8,9])

    # filter for rep_names
    rep_names = ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]
    logging.info(f"Filtering for rep_names: {rep_names}")
    df0 =  df0[df0['repeat'].isin(rep_names)]

    # save to new dataframe
    df1 = pd.DataFrame()
    df1['chr'] = df0['chr']
    # set start positions depending on strand
    df1['pos'] = df0.apply(lambda x: x['end'] if x['strand'] != '+' else x['start'], axis=1)
    df1['l1_name'] = df0['repeat']

    return df1

def main():

    # setup logging
    logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

    df = read_rmsk()
    l1_pos = set()

    for (_, chrom, pos) in df[['chr', 'pos']].itertuples():
        l1_pos.update([Interval(chrom, pos, pos)])
        
    genome = Genome(snakemake.input['genome'])
    xx = list(ig.windows_overlapping_intervals(genome, l1_pos, 750, 250))

    refl1 = pd.DataFrame.from_records((x.as_tuple() for x in xx), columns=['chr', 'start', 'end']).set_index(['chr', 'start', 'end'])
    refl1['reference_l1hs_l1pa2_6'] = True

    refl1.to_csv(snakemake.output[0])

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