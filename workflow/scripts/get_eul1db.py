#!/usr/bin/env python
# adapted from Ricardo's jupyter notebook from rf pipeline: pipeline/eul1db/to_windows/_h/main.ipynb
__author__ = 'Michael Cuoco'

# TODO: add comments 

import pandas as pd
from pyslavseq.genome import Interval, Genome
import pyslavseq.genome.interval_generator as ig
    
def main():
    url = "http://eul1db.unice.fr/UserLists/DATA/downloads/SRIP.txt"
    df0 = pd.read_csv(url, sep="\t", skiprows=5, dtype={'chromosome':str})
    df0['chromosome'] = 'chr' + df0['chromosome']

    df = df0[(df0['lineage']=='germline') & (df0['study_id'].isin(['Ewing2010', 'Ewing2011', 'Stewart2011', 'Beck2010', 'Brouha2002', 'Iskow2010'])) & (df0['g_start']==df0['g_stop'])]

    bed = df[['chromosome', 'g_start', 'g_stop']]
    bed = bed.rename(columns={'chromosome': 'chr', 'g_start': 'start', 'g_stop': 'stop'})
    bed['start'] -= 1

    bed.to_csv(snakemake.output[0], sep="\t", header=True)

if __name__ == '__main__':
    main()