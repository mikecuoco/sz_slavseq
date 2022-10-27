#!/usr/bin/env python
__author__ = 'Apua Paquola, Michael Cuoco'

import argparse
import pickle
import gzip
import pandas as pd
from pyslavseq.genome import Genome
from get_windows import make_genome_windows
import sys, gc, traceback

def flank_features(df):
    for i in range(1, 8):
        flank_size = 2**i
        field_name = 'flank_' + str(flank_size) + '_max_reads'
        df[field_name] = df['all_reads.count'].rolling(window=2*flank_size + 1, center=True, min_periods=1).max().fillna(0)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bgz', required=True)
    parser.add_argument('--chromsizes', required=True)
    parser.add_argument('--outfile', required=True)
    args=parser.parse_args()

    emptydf = make_genome_windows(args.chromsizes)

    df = pd.read_csv(args.bgz, compression='gzip', sep='\t', index_col=[0, 1, 2])
    arcdf = df[['all_reads.count']]
    arcdf = pd.merge(emptydf, arcdf, how='outer', left_index=True, right_index=True).fillna(0)
    flank_features(arcdf)
    arcdf = arcdf[arcdf['all_reads.count'] > 0]
    arcdf = arcdf.drop('all_reads.count', axis=1)

    df = pd.merge(df, arcdf, how='inner', left_index=True, right_index=True)

    f = gzip.open(args.outfile,'wb')
    pickle.dump(df, f)
    f.close()

if __name__ == '__main__':
    main()
