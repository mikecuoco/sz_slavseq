#!/usr/bin/env python
# adapted from Ricardo's jupyter notebooks from rf pipeline's somaticMutationSummary step
__author__ = "Rohini Gadde"

import itertools
import numpy as np
import pandas as pd
import os
import re
import pyslavseq
from pyslavseq.genome import Genome, Interval, overlaps

from pathlib import Path

from IPython.display import clear_output

import warnings
warnings.filterwarnings("ignore")

# Merge predictions from random forest
def get_prediction():
   
    print("\tMerging prediction of Fold-%d ..." % k)
    df_open =  pd.read_csv(fileFold + 'Testing_y_pred.csv', sep=";")   
    df_open.set_index( ['chrom','start','end'] ,inplace=True) 
    
    return df_open

def sample_id(s):
    new = s.str.split(pat = "_", n = 1, expand = True) 

    return new[0], new[1]

def get_summary(sample,bulk=False):        
    
    df.sort_index(ascending=True, inplace=True)
        
    df.reset_index(inplace=True)

    if bulk: 
        sample_name, cell = sample.split("_",1)[0],'bulk'
    else:
        sample_name, cell = sample_id (df['cell_id']) 
        
    df['end'] = sample_name
    df['cell_id'] = cell

    df.rename({'end': 'sample'}, axis=1, inplace=True)
 
    return df

def main():
    inDirs = [Path(f) for f in snakemake.input if re.search(f'fold_{k}', f)]

    header=['chrom','start','sample','cell_id','Y','Y_pred','all_reads_count','KNRGL_proba','OTHER_proba','RL1_proba']
    
    for sample in inDirs:

        df = pd.DataFrame(data=None)
        
        print("\n")
        print("Processing Sample ID ({}) ... \n" .format( sample ) )
            
        for k in range(0,4):
            
            fileFold = sample + "/"

            df = df.append( get_prediction() )
        
        if 'bulk' in sample: bulk = True 
        else: bulk = False
        
        df0 = get_summary(sample,bulk)
            
        df0[header].to_csv(sample + '.csv', sep="\t",index=False,header=True )

if __name__ == '__main__':
    main()