#!/jhpce/shared/jhpce/core/python/3.8.3/bin/python3
__author__ = 'Ricardo S Jacomini, Michael Cuoco'

import argparse
import sys,gc
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

import pandas as pd
import pickle
import functools
import re
import numpy as np

import pyslavseq
from pyslavseq.genome import Interval, Genome, is_within
import pyslavseq.genome.interval_generator as ig

def interval_hash(iv, window_size):
    return iv.chrom + "_" + str(iv.start // window_size)
    
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

@functools.lru_cache()
def get_reference_l1():
    df = pd.read_csv(snakemake.input.ref_l1, index_col=[0,1,2])
    return df

@functools.lru_cache()
def get_eul1db():
    df = pd.read_csv(snakemake.input.non_ref_l1, index_col=[0,1,2])
    return df

# @functools.lru_cache(maxsize=1500)
def get_cell_features(cell_id):
    fn = "../../_m/%s.pickle.gz" % cell_id
    df = pd.read_pickle(fn).merge(get_eul1db(), left_index=True, right_index=True, how="left")\
    .merge(get_reference_l1(), left_index=True, right_index=True, how="left")\
    .fillna({'in_eul1db':False, 'reference_l1hs_l1pa2_6':False})
    df['cell_id'] = cell_id
    return df

def cells_from_sample(sample_files, dna_type='mda'):
    for fn in sample_files:
        if dna_type == 'bulk': 
            if 'bulk' not in fn: continue
            m = re.search('/([^/]+?)\.pickle\.gz$', fn)
            yield m.group(1)
        else:
            if 'bulk' in fn: continue
            m = re.search('/([^/]+?)\.pickle\.gz$', fn)
            yield m.group(1)

def my_fold(df, window_size, num_folds):
    genome = Genome(snakemake.input.chromsizes)
    d = interval_dict(genome, window_size)

    for x in df.index:
        # index is ['chrom', 'start', 'end', 'cell_id']
        iv = Interval(x[0], x[1], x[2])
        yield interval_fold_assignment(d, iv, window_size, num_folds)

def my_label(df):
    for x in df[['in_eul1db','reference_l1hs_l1pa2_6']].itertuples():
        if x.reference_l1hs_l1pa2_6:
            yield 'RL1'
        elif x.in_eul1db:
            yield 'KNRGL'
        else:
            yield 'OTHER'
            
def get_train_and_test_data(df, fold, min_reads):

    features = [x for x in df.columns if not any(x.endswith(y) for y in [
        'chrom', 'start', 'end', 'cell_id', 'in_eul1db', 'reference_l1hs_l1pa2_6',
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


def folds (sample, dna_type):
 
    print(sample)
 
    if dna_type == 'bulk':
        x = [x for x in cells_from_sample(sample,dna_type)]
        df = get_cell_features(x[0])
    else:
        df = pd.concat([get_cell_features(x).reset_index() for x in cells_from_sample(sample,dna_type)]).sort_values(['chrom', 'start', 'end', 'cell_id']).set_index(['chrom', 'start', 'end', 'cell_id'])
        
    df['fold'] = list(my_fold(df, window_size, num_folds))

    for fold in range(num_folds):
        out_dir = '%s_%s/fold_%d' % (sample, dna_type, fold)
        
        print(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        
        X_train, Y_train, X_test, Y_test = get_train_and_test_data(df, fold, min_reads)
        
        X_train.to_pickle(out_dir + '/X_train.pickle.gz', compression='gzip')
        Y_train.to_pickle(out_dir + '/Y_train.pickle.gz', compression='gzip')
        X_test.to_pickle(out_dir + '/X_test.pickle.gz', compression='gzip')
        Y_test.to_pickle(out_dir + '/Y_test.pickle.gz', compression='gzip')
        
        #X_train.to_csv(out_dir + '/X_train.csv')

def print_err(msg):
    sys.stderr.write(msg + '\n')

def main():
    
    global window_size, num_folds, min_reads
    window_size = snakemake.params.window_size
    num_folds = snakemake.params.num_folds
    min_reads = snakemake.params.min_reads

    sample_files = snakemake.input.samples

    df = pd.concat([get_cell_features(x).reset_index() for x in cells_from_sample(sample)]).sort_values(['chrom', 'start', 'end', 'cell_id']).set_index(['chrom', 'start', 'end', 'cell_id'])

    df['fold'] = list(my_fold(df, window_size, num_folds))

    for fold in range(num_folds):
        out_dir = 'fold_%d' % fold
        print(out_dir)
        os.makedirs(out_dir, exist_ok=True)
    
        X_train, Y_train, X_test, Y_test = get_train_and_test_data(df, fold, min_reads)
    
        X_train.to_pickle(out_dir + '/X_train.pickle.gz', compression='gzip')
        X_test.to_pickle(out_dir + '/X_test.pickle.gz', compression='gzip')
        Y_train.to_pickle(out_dir + '/Y_train.pickle.gz', compression='gzip')
        Y_test.to_pickle(out_dir + '/Y_test.pickle.gz', compression='gzip')
    
    #X_train.to_csv(out_dir + '/X_train.csv')
    #Y_train.to_csv(out_dir + '/Y_train.csv', header=True)
    #X_test.to_csv(out_dir + '/X_test.csv')
    #Y_test.to_csv(out_dir + '/Y_test.csv', header=True)
        
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=False)
    parser.add_argument('--dna_type', required=True)
    
    args=parser.parse_args()
    
    try:
        folds(args)
    except: # catch *all* exceptions
        print_err( "Message : %s" % sys.exc_info()[0])
        print_err( "%s" % sys.exc_info()[1])
        print_err( "%s" % sys.exc_info()[2])
    finally:
        # cleanup code in here
        gc.collect()
