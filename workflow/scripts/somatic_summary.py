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

def get_summary(donor, dna_type, df, bulk=False):        
    
    df.sort_values(by='start', ascending=True, inplace=True)
        
    # df.reset_index(inplace=True)

    if bulk: 
        sample_name, cell = donor, "bulk"
    else:
        sample_name, cell = donor, dna_type
        
    df['end'] = sample_name
    df['cell_id'] = cell

    df.rename({'end': 'sample'}, axis=1, inplace=True)
 
    return df

def main():
    # input is in format "results/train_test/{{donor}}/{{dna_type}}/{fold}"
    
    # get input directory containing each fold directory
    inDirs = [Path(d).parent.resolve() for d in snakemake.input]
    inDirs = np.unique(np.array(inDirs)) # get only unique directories

    dna_types = [d.name for d in inDirs]
    donors = [Path(d).parent.name for d in inDirs]

    outDirs = [str(d).replace("train_test", "somatic_summary") for d in inDirs]
    
    for i in range(0, len(inDirs)):
        sample = donors[i] + "_" + dna_types[i]
        inDir = str(inDirs[i])

        df = pd.DataFrame(data=None)
        
        print("\n")
        print("Processing Sample ID ({}) ... \n" .format(sample))
            
        for k in range(0, snakemake.params.num_folds):
            
            pred = inDir + "/" + "fold_" + str(k) + "/" + "Testing_y_pred.csv"
            
            if not os.path.exists(pred): # skip folds with no equal classes
                continue
            
            print("\tMerging prediction of Fold-%d ..." % k)
            df_open = pd.read_csv(pred)
            df = df.append(df_open)
        
        if 'bulk' in sample: bulk = True
        else: bulk = False
        
        df0 = get_summary(donors[i], dna_types[i], df, bulk)
        
        if not os.path.exists(outDirs[i]):
            os.mkdir(outDirs[i])
            
        df0.to_csv(outDirs[i] + "/" + "Merged_y_pred.csv", index=False, header=True)

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
