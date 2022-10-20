#!/usr/bin/env python
# adapted from Ricardo's jupyter notebook from rf pipeline: pipeline/eul1db/to_windows/_h/main.ipynb
__author__ = 'Michael Cuoco, Rohini Gadde'

import pandas as pd
from pyslavseq.genome import Interval, Genome
import pyslavseq.genome.interval_generator as ig
import snakemake as sm
import logging
import sys, gc, traceback

def read_rmsk(rmsk_outfile=snakemake.input['rmsk']):
    """
    Read the repeatmasker output table and return locations of L1HS and L1PA2-6
    """
    # read the rmsk file
    df0 = pd.read_csv(rmsk_outfile, skiprows=3, delim_whitespace=True, 
        names=["chrom", "start", "end", "strand", "repeat"], usecols=[4,5,6,8,9])

    # filter for rep_names
    rep_names = ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]
    logging.info(f"Filtering for rep_names: {rep_names}")
    df0 =  df0[df0['repeat'].isin(rep_names)]

    # save to new dataframe
    df1 = pd.DataFrame()
    df1['chrom'] = df0['chrom']
    # set start positions depending on strand
    df1['pos'] = df0.apply(lambda x: x['end'] if x['strand'] != '+' else x['start'], axis=1)
    df1['l1_name'] = df0['repeat']

    return df1

def make_intervals(df, chr_sizes = snakemake.input['chromsizes'], outfile=snakemake.output[0]):
    """
    Make L1 positions into genomic intervals
    """
    # TODO: add interval size and overlap as parameters
    # TODO: move this function out of this file, to be accessible to other scripts
    l1_pos = set()

    for (_, chrom, pos) in df[['chrom', 'pos']].itertuples():
        l1_pos.update([Interval(chrom, pos, pos)])
        
    genome = Genome(chr_sizes)
    xx = list(ig.windows_overlapping_intervals(genome, l1_pos, 750, 250))

    refl1 = pd.DataFrame.from_records((x.as_tuple() for x in xx), columns=['chrom', 'start', 'end']).set_index(['chrom', 'start', 'end'])
    refl1['reference_l1hs_l1pa2_6'] = True

    refl1.to_csv(outfile)

if __name__ == '__main__':
    try:
        df = read_rmsk()
        make_intervals(df)

    except:  # catch *all* exceptions
        sys.stderr = open(snakemake.log[0], 'w')
        traceback.print_exc()
        sys.stderr.close()

    finally:
        # cleanup code in here
        gc.collect()