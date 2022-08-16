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
import subprocess

import pyslavseq
from pyslavseq.genome import Genome, Interval, overlaps

from pathlib import Path

# from IPython.display import clear_output

import warnings
warnings.filterwarnings("ignore")

def get_summary(donor, dna_type, df):        
    
    df.sort_values(by='start', ascending=True, inplace=True)
        
    # df.reset_index(inplace=True)

    sample_name, cell = donor, dna_type
        
    df['end'] = sample_name
    df['cell_id'] = cell

    df.rename({'end': 'sample'}, axis=1, inplace=True)
 
    return df

class Window(Interval):
    """A Window with position (a string, start is 0-based, end is 1-based) 
    where of intervals come from random forest and slav-seq features
    """
    def __init__(self, chrom=None, start=None, end=None,proba=None,reads=None):
        assert(start <= end)
        super().__init__(chrom, start, end)
        self.proba=proba
        self.reads=reads

        self.interval = pd.Interval(start,end,closed='both')
        
    def __str__(self):
        return str(super().__str__()) + "\t Score: " + str(self.proba) + '\t Reads: ' + str(self.reads)
    
    def as_tuple(self):
        return self.chrom, self.start, self.end, self.interval, self.proba, self.reads
    
class Cell():
    """A collection of Windows
    """
    def __init__(self, cell=None, listWindows=None):
        self.cell=cell

        if listWindows is None:
            listWindows = []
        self.windows = listWindows
        self.windows_max= []

    def __str__(self):
        return self.cell
    
    def append(self,w):
        self.windows.append(w)
    
    def __printlist(self,df):
        print( self.cell + "\n")
        
        if len(df) == 0: print("There is no window !")
        else:
            if len(df) == 1:  print(df)
            else:
                for w in sorted(df):
                    print(w)
            
    def print_windows(self):
        self.__printlist(self.windows)
    
    def print_windows_max(self):
        self.__printlist(self.windows_max)

    def to_dataframe(self,sample=None):
        
        if len(self.windows_max) > 1:
            df = pd.DataFrame(data=None)

            for w in sorted(self.windows_max):
                df = df.append( [ w.as_tuple() ] )
        else:
            import numpy as np
            df = np.transpose((pd.DataFrame(data=w1.as_tuple())))

        df.columns = ['chrom','start','end','intv', 'proba_max', 'reads_max']
       
        df['cell_id'] = self.cell
        df['sample'] = sample
        
        return df.reset_index(drop=True)
        
    def windows_overlap(self):
        w = sorted(self.windows)
        
        i = 1
        
        if len(w) > 1:
            self.windows_max.append( w[0] )

            while i <= (len(w)-1):

                x = is_overlap( self.windows_max[-1], w[i])

                if x != w[i]: 
                    self.windows_max.pop()

                self.windows_max.append(x)

                i += 1
        else:
            self.windows_max.append(w)
                    
def is_overlap(iv1, iv2):
    """Returns new Window with overlapping part of two intervals if they overlap.
    """
    if overlaps_local(iv1, iv2):
        return Window(iv1.chrom, 
                      min(iv1.start, iv2.start), 
                      max(iv1.end, iv2.end),
                      max(iv1.proba, iv2.proba),
                      max(iv1.reads, iv2.reads))    
    else:
        return iv2
    
def overlaps_local(iv1, iv2):
    """Returns a boolean indicating whether intervals overlap.
    """
    return iv1.chrom == iv2.chrom and iv1.interval.overlaps(iv2.interval)

def my_cluster_id(df):
    
    prev_iv = None
    cluster_id = -1
    for x in df[['chrom','start','end']].itertuples():
        iv = Interval(x[1], x[2], x[3])
        if prev_iv is None or not overlaps(iv, prev_iv):
            cluster_id += 1
        yield cluster_id
        prev_iv = iv
        
def my_cluster_id_cell(df):
    prev_iv = None
    prev_cell = None
    cluster_id = -1
    for x in df[['chrom','start','end', 'sample', 'cell_id']].itertuples():
        iv = Interval(x[1], x[2], x[3])
        if prev_iv is None or not overlaps(iv, prev_iv) \
        or prev_cell is None or prev_cell != (x[4], x[5]):
            cluster_id += 1
        yield cluster_id
        prev_iv = iv
        prev_cell = (x[4], x[5]) 

def get_precision_position(fileFold, metric):
    i = -1; a = ""
    while metric not in a:
        i += 1
        a = open(fileFold, "r").readlines()[i]
    return i

def get_results_metrics(inDirs, samples):
    
    header = ['sample','fold','train_precision','train_recall','train_f1','test_precision','test_recall','test_f1']

    a  = np.zeros(shape=(1,8))
    
    df = pd.DataFrame(data=None,columns=header)
    X2 = pd.DataFrame(data=a,columns=header)

    for i in range(0, len(inDirs)):
        tmp = X2.copy()

        tmp['sample'] = samples[i]
        print("Loading dataset: " + samples[i])

        for k in range(0, snakemake.params.num_folds):

            tmp['fold'] = k
            fn = "/Training_report_performance.txt"
            tFile ='train'
        
            fileFold =  str(inDirs[i]) + "/fold_" + str(k) + "/" + fn

            precision = get_precision_position(fileFold, "precision")
            recall = get_precision_position(fileFold, "recall")
            f1 = get_precision_position(fileFold, "f1")

            lines = [precision, recall, f1]

            for line in lines:
                a = open(fileFold, "r").readlines()[line]
                a = a.replace("\n", "")

                if line == lines[0]: tmp[tFile + '_precision'] = float(a.replace("precision ", ""))
                if line == lines[1]: tmp[tFile + '_recall'] = float(a.replace("recall ", ""))
                if line == lines[2]: tmp[tFile + '_f1'] = float(a.replace("f1 ", ""))

            fn = "/Testing_report_performance.txt"
            tFile ='test'
                
            fileFold = str(inDirs[i]) + "/fold_" + str(k) + "/" + fn

            for line in lines:
                a = open(fileFold, "r").readlines()[line]
                a = a.replace("\n", "")

                if line == lines[0]: tmp[tFile + '_precision'] = float(a.replace("precision ", ""))
                if line == lines[1]: tmp[tFile + '_recall'] = float(a.replace("recall ", ""))
                if line == lines[2]: tmp[tFile + '_f1'] = float(a.replace("f1 ", ""))
            
            df = df.append(tmp)

    return df

def main():
    # input is in format "results/train_test/{{donor}}/{{dna_type}}/{fold}"
    # output is in "results/somatic_summary/{{donor}}/{{dna_type}}/" directory

    # get input directory containing each fold directory
    inDirs = [Path(d).parent.resolve() for d in snakemake.input]
    inDirs = np.unique(np.array(inDirs)) # get only unique directories

    dna_types = [d.name for d in inDirs]
    donors = [Path(d).parent.name for d in inDirs]
    samples = [donors[i] + "_" + dna_types[i] for i in range(0, len(inDirs))]

    outDirs = [str(d).replace("train_test", "somatic_summary") for d in inDirs]
    
    window_size = snakemake.params.window_size
    prob = snakemake.params.min_prob
    reads_count = snakemake.params.min_reads

    slavseq_sz = pd.DataFrame(data=None)

    for i in range(0, len(inDirs)):
        sample = samples[i]
        inDir = str(inDirs[i])

        df_merge = pd.DataFrame(data=None)
        
        print("\n")
        print("Processing Sample ID ({}) ... \n" .format(sample))
            
        for k in range(0, snakemake.params.num_folds):
            
            pred = inDir + "/" + "fold_" + str(k) + "/Testing_y_pred.csv"
            
            if not os.path.exists(pred): # skip folds with no equal classes
                continue
            
            # Merge predictions from each fold into one
            print("\tMerging prediction of Fold-%d ..." % k)
            df_open = pd.read_csv(pred)
            df_merge = df_merge.append(df_open)
        
        # if 'bulk' in sample: bulk = True
        # else: bulk = False
        
        # add sample and cell labels
        df_label = get_summary(donors[i], dna_types[i], df_merge)
        
        if not os.path.exists(outDirs[i]):
            os.mkdir(outDirs[i])

        df_label.to_csv(outDirs[i] + "/Merged_y_pred.csv", index=False, header=True)

        
        df_sz = pd.DataFrame(data=None)
        
        # filter for potential somatic insertions
        cond2 = df_label['KNRGL_proba'] >= prob
        cond3 = df_label['Y'] == 'OTHER' 
        cond4 = df_label['Y_pred'] == 'KNRGL'
        cond5 = df_label['all_reads_count'] >= reads_count

        df_filter = df_label[ cond2 & cond3 & cond4 & cond5 ].reset_index(drop=True)
        
        # if bulk: df_filter['cell_id'] = 'bulk'
        
        cells = set(df_filter['cell_id'])
        
        if len(df_filter) > 0:
            for c in cells:

                print(sample, "-", c)

                c1 = Cell(cell=c)
                
                # drop rows with different cell type
                cond1 = df_filter['cell_id'] == c 
                df_cell = df_filter[ cond1 ].reset_index(drop=True)

                for r in range(0, len(df_cell)):
                    w = df_cell.iloc[r]

                    w1 = Window(w.chrom,w.start,w.start+window_size,w.OTHER_proba,w.all_reads_count)

                    c1.append(w1)
                c1.windows_overlap() # return overlapping windows

                df_sz = df_sz.append(c1.to_dataframe(sample), ignore_index=True)         
                
            slavseq_sz = slavseq_sz.append(df_sz, ignore_index=True)
    
            df = slavseq_sz.sort_values(['sample', 'cell_id', 'chrom', 'start', 'end']).copy()
            
            # type conversions
            df['start'] = df.start.astype('int32') 
            df['end'] = df.end.astype('int32') 
            df['window_id'] = list(my_cluster_id_cell(df))

            df = df.sort_values(['window_id'], ascending=True).groupby('window_id').first()

            df.drop(['intv'], axis=1, inplace=True) # drop interval [start, end]

            df = df.sort_values(['chrom', 'start', 'end'])
            df['cluster_id'] = list(my_cluster_id(df)) 
            
            cs = df.groupby('cluster_id').size().to_frame().rename(columns={0:'number_cells_detected'})
            df = df.merge(cs, left_on='cluster_id', right_index=True) # add new col # cells detected

            # filter df for given columns
            df = df[
                ['chrom', 'start', 'end', 'sample', 'cell_id', 'proba_max',
                'reads_max', 'number_cells_detected','cluster_id']
                ]
            
            df = df.rename(columns={'proba_max':'confidence_score', 'reads_max': 'number_reads'})
            df.to_csv(outDirs[i] + "/slavseq_sz-intersections-cluster.csv", index=True, header=True)

            df = df[
                (df['number_cells_detected'] <= 5) ].sort_values(
                    ['sample','cell_id','confidence_score'], 
                    ascending=[False, False, False])

            df.to_csv(outDirs[i] + "/somatic_candidates-cluster.csv", index=True, header=True)
        
        else:
            # output empty csvs with only header if no insertions found
            slavseq_sz = pd.DataFrame(columns = ['chrom', 'start', 'end', 'sample', 'cell_id', 'confidence_score',
                'number_reads', 'number_cells_detected','cluster_id'])
            slavseq_sz.to_csv(outDirs[i] + "/slavseq_sz-intersections-cluster.csv", header=True)
            slavseq_sz.to_csv(outDirs[i] + "/somatic_candidates-cluster.csv", header=True)
    
    
    df_metrics = get_results_metrics(inDirs, samples)
    df_metrics.to_csv(outDirs[i] + "/Cross_validation_metrics.csv", index=False, header=True)


    # Merge results of "Testing_y_pred.csv" for ALL folds for ALL samples
    subprocess.run("head -1 " + outDirs[0] + "/Merged_y_pred.csv" + "  >" + outDirs[0] + 
        '/slavseq_sz_no_filter.csv' + " 2>"  + outDirs[0] + '/slavseq_sz_no_filter.stderr', 
        shell=True)

    for out in outDirs:
        subprocess.run("tail -n +2 " + out + "/Merged_y_pred.csv" + " >>" + out + 
        '/slavseq_sz_no_filter.csv' + " 2>>" + out + '/slavseq_sz_no_filter.stderr', 
            shell=True)

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
