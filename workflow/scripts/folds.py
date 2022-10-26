#!/usr/bin/env python
__author__ = 'Ricardo S Jacomini, Michael Cuoco'

import sys, gc, traceback
import os
os.environ['OPENBLAS_NUM_THREADS'] = str(snakemake.threads)

import pandas as pd
import functools
import re
import numpy as np
from pathlib import Path
from pyslavseq.genome import Interval, Genome, is_within
import pyslavseq.genome.interval_generator as ig

def interval_hash(iv, window_size):
    return str(iv.chrom) + "_" + str(iv.start // window_size)
    
@functools.lru_cache()
def interval_dict(genome, window_size):
    d = dict()
    for ii, iv in enumerate(ig.windows_in_genome(genome, window_size, window_size)):
        d[interval_hash(iv, window_size)] = (ii, iv)

    return d

def interval_fold_assignment(d, iv, window_size, num_folds):
    # Assigns intervals to folds (for cross validation)
    # Assign None to intervals that cross fold boundaries
    (interval_index, target_interval) = d[interval_hash(iv, window_size)]
    if is_within(iv, target_interval):
        return interval_index % num_folds
    else:
        return -1

def read_rmsk(rmsk_outfile=snakemake.input.ref_l1):
    """
    Read the repeatmasker output table and return locations of L1HS and L1PA2-6
    """
    # read the rmsk file
    df0 = pd.read_csv(rmsk_outfile, skiprows=3, delim_whitespace=True, 
        names=["chrom", "start", "end", "strand", "repeat"], usecols=[4,5,6,8,9])

    # filter for rep_names
    rep_names = ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]
    # logging.info(f"Filtering for rep_names: {rep_names}")
    df0 =  df0[df0['repeat'].isin(rep_names)]

    # save to new dataframe
    df1 = pd.DataFrame()
    df1['chrom'] = df0['chrom']
    # set start positions depending on strand
    df1['start'] = df0.apply(lambda x: x['end'] if x['strand'] != '+' else x['start'], axis=1)
    df1['end'] = df1['start']
    df1['start'] -= 1 # make zero-based

    return df1

def make_windows(df, field):
    l1_pos = set()

    for (_, chrom, start, end) in df[['chrom', 'start', 'end']].itertuples():
        l1_pos.update([Interval(chrom, start, end)])

    genome = Genome(snakemake.input.chromsizes)
    xx = list(ig.windows_overlapping_intervals(genome, l1_pos, 750, 250))
    
    l1_df = pd.DataFrame.from_records((x.as_tuple() for x in xx), columns=['chrom', 'start', 'end']).set_index(['chrom', 'start', 'end'])
    l1_df[field] = True 
    
    return l1_df

@functools.lru_cache()
def get_reference_l1():
    rmsk_df = read_rmsk()
    l1_df = make_windows(rmsk_df, "reference_l1hs_l1pa2_6")
    # df = pd.read_csv(snakemake.input.ref_l1[0], index_col=[0,1,2])
    return l1_df

@functools.lru_cache()
def get_non_ref_db():
    nr_df = pd.read_csv(snakemake.input.non_ref_l1, sep="\t", names=["chrom", "start", "end"])
    l1_df = make_windows(nr_df, "in_NRdb")
    # df = pd.read_csv(snakemake.input.non_ref_l1[0], index_col=[0,1,2])
    return l1_df

# @functools.lru_cache(maxsize=1500)
def get_cell_features(fn, cell_id):
    # merge reference and non-reference germline L1s
    df = pd.read_pickle(fn).merge(get_non_ref_db(), left_index=True, right_index=True, how="left")\
        .merge(get_reference_l1(), left_index=True, right_index=True, how="left")\
        .fillna({'in_NRdb':False, 'reference_l1hs_l1pa2_6':False})
        
    df['cell_id'] = cell_id
    
    return df

def cells_from_sample(sample_files):
    for fn in sample_files:
        m = re.search('/([^/]+?)\.pickle\.gz$', fn)
        yield fn, m.group(1)

def my_fold(df, window_size, num_folds):
    genome = Genome(snakemake.input.chromsizes)
    d = interval_dict(genome, window_size)

    for x in df.index:
        # index is ['chrom', 'start', 'end', 'cell_id']
        iv = Interval(x[0], x[1], x[2])
        yield interval_fold_assignment(d, iv, window_size, num_folds)

def my_label(df):
    for x in df[['in_NRdb','reference_l1hs_l1pa2_6']].itertuples():
        if x.reference_l1hs_l1pa2_6:
            yield 'RL1'
        elif x.in_NRdb:
            yield 'KNRGL'
        else:
            yield 'OTHER'
            
def get_train_and_test_data(df, fold, min_reads):

    features = [x for x in df.columns if not any(x.endswith(y) for y in [
        'chrom', 'start', 'end', 'cell_id', 'in_NRdb', 'reference_l1hs_l1pa2_6',
        'fold', '_peak_position', '_en_motif', '_te_strand'])]
        
    train_examples = (df['all_reads.count']>=min_reads) & (df['fold']>=0) & (df['fold']!=fold)
    test_examples = (df['all_reads.count']>=min_reads) & (df['fold']==fold)
    
    # train features
    X_train = df[train_examples][features].fillna(0)
    for field in ['all_reads.median_r1r2_distance', 'yayg_reads.median_r1r2_distance', 
                  'yg_reads.median_r1r2_distance', 'all_reads.secondary.median_distance', 
                  'yayg_reads.secondary.median_distance', 'yg_reads.secondary.median_distance']:
        X_train[field] = np.minimum(X_train[field], 4e9)
    
    # train labels
    Y_train = pd.Series(my_label(df[train_examples]), index=df[train_examples].index)
    
    # test features
    X_test = df[test_examples][features].fillna(0)
    for field in ['all_reads.median_r1r2_distance', 'yayg_reads.median_r1r2_distance', 
                  'yg_reads.median_r1r2_distance', 'all_reads.secondary.median_distance', 
                  'yayg_reads.secondary.median_distance', 'yg_reads.secondary.median_distance']:
        X_test[field] = np.minimum(X_test[field], 4e9)
    
    # test labels
    Y_test = pd.Series(my_label(df[test_examples]), index=df[test_examples].index)
    
    return X_train, Y_train, X_test, Y_test


def folds(sample_files):
 
    if dna_type == 'bulk':
        x = [(fn, cell_id) for fn, cell_id in cells_from_sample(sample_files)]
        df = get_cell_features(*x[0])
    else:
        df = pd.concat([get_cell_features(fn, cell_id).reset_index() for fn, cell_id in cells_from_sample(sample_files)]).sort_values(['chrom', 'start', 'end', 'cell_id']).set_index(['chrom', 'start', 'end', 'cell_id'])

    df['fold'] = list(my_fold(df, fold_window, num_folds))

    for fold in range(num_folds):
        # get directory for output files of this fold
        fold_dir = set([str(Path(f).parent) for f in snakemake.output if re.search(f'fold_{fold}', f)]).pop()

        X_train, Y_train, X_test, Y_test = get_train_and_test_data(df, fold, min_reads)
        
        X_train.to_pickle(fold_dir + '/X_train.pickle.gz', compression='gzip')
        Y_train.to_pickle(fold_dir + '/Y_train.pickle.gz', compression='gzip')
        X_test.to_pickle(fold_dir + '/X_test.pickle.gz', compression='gzip')
        Y_test.to_pickle(fold_dir + '/Y_test.pickle.gz', compression='gzip')

def print_err(msg):
    sys.stderr.write(msg + '\n')
        
if __name__ == '__main__':

    global fold_window, num_folds, min_reads, dna_type
    fold_window = snakemake.params.fold_window
    num_folds = snakemake.params.num_folds
    min_reads = snakemake.params.min_reads
    dna_type = snakemake.wildcards.dna_type

    sample_files = snakemake.input.samples

    # not sure if try/except is necessary
    try:
        folds(sample_files)
    except: # catch *all* exceptions
        sys.stderr = open(snakemake.log[0], 'w')
        traceback.print_exc()
        sys.stderr.close()
    finally:
        # cleanup code in here
        gc.collect()

