#!/usr/bin/env python
# adapted from Ricardo's jupyter notebooks from rf pipeline's somaticMutationSummary step
__author__ = "Rohini Gadde"

import sys,gc
import traceback
import os
import itertools
import numpy as np
import pandas as pd
import re
import pyslavseq
from pyslavseq.genome import Genome, Interval, overlaps

from pathlib import Path

# from IPython.display import clear_output

import warnings
warnings.filterwarnings("ignore")

def get_prediction(dirPath, header):
   
    # Add this line to handle folds with no equal classes
    if not os.path.exists(dirPath + 'Testing_y_pred.csv'):
        df_open = pd.DataFrame(columns = header)
        # df_open.set_index( ['chrom','start','end'], inplace=True)
    else:
        print("\tMerging prediction of Fold-%d ..." % k)
        df_open =  pd.read_csv(dirPath + 'Testing_y_pred.csv', sep=";")
        df_open.set_index( ['chrom','start','end'], inplace=True)
    
    return df_open

# def sample_id(s):
#     new = s.str.split(pat = "_", n = 1, expand = True) 

#     return new[0], new[1]

def get_summary(sample, df, bulk=False):        
    
    df.sort_index(ascending=True, inplace=True)
        
    df.reset_index(inplace=True)

    if bulk: 
        sample_name, cell = sample, 'bulk'
    else:
        sample_name, cell = sample, df['cell_id'] 
        
    df['end'] = sample_name
    df['cell_id'] = cell

    df.rename({'end': 'sample'}, axis=1, inplace=True)
 
    return df

def main():
    inDirs = [Path(d).name for d in snakemake.input if re.search(f'fold_[0-9]', d)]

    header=['chrom','start','sample','cell_id','Y','Y_pred','all_reads_count','KNRGL_proba','OTHER_proba','RL1_proba']
    
    for sample in inDirs:

        df = pd.DataFrame(data=None)
        
        print("\n")
        print("Processing Sample ID ({}) ... \n" .format( sample ) )
            
        for k in range(0, snakemake.params.num_folds):
            
            fileFold = sample + "/"
            
            df = df.append( get_prediction(fileFold, header) )
        
        if 'bulk' in sample: bulk = True 
        else: bulk = False
        
        df0 = get_summary(sample, df, bulk)
        # TODO: change output path
        df0[header].to_csv(sample + '.csv', sep="\t", index=False, header=True )

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
