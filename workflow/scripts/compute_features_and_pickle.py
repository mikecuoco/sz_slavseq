#!/usr/bin/env python
__author__ = 'Apua Paquola'

#import os
#os.environ['OPENBLAS_NUM_THREADS'] = '1'

import argparse
import pickle
import gzip

import pandas as pd
from pyslavseq.genome import interval_generator as ig
from pyslavseq.genome import Genome


def genome_empty_df(chromsizes):
    genome = Genome(chromsizes)
    wg_iter = ig.windows_in_genome(genome, 750, 250)
    df = pd.DataFrame.from_records((xx.as_tuple() for xx in wg_iter), columns=['chrom', 'start', 'end'])
    return df.set_index(['chrom', 'start', 'end'])


def flank_features(df):
    for i in range(1, 8):
        flank_size = 2**i
        field_name = 'flank_' + str(flank_size) + '_max_reads'
        df[field_name] = df['all_reads.count'].rolling(window=2*flank_size + 1, center=True, min_periods=1).max().fillna(0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--chromsizes', required=True)
    parser.add_argument('--output', required=True)
    args=parser.parse_args()

    emptydf = genome_empty_df(args.chromsizes)

    df = pd.read_csv(args.input, compression='gzip', sep='\t', index_col=[0, 1, 2])
    arcdf = df[['all_reads.count']]
    arcdf = pd.merge(emptydf, arcdf, how='outer', left_index=True, right_index=True).fillna(0)
    flank_features(arcdf)
    arcdf = arcdf[arcdf['all_reads.count'] > 0]
    arcdf = arcdf.drop('all_reads.count', axis=1)

    df = pd.merge(df, arcdf, how='inner', left_index=True, right_index=True)

    f = gzip.open(args.output,'wb')
    pickle.dump(df, f)
    f.close()


if __name__ == '__main__':
    main()
