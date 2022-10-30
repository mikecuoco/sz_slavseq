#!/usr/bin/env python
# adapted from Ricardo's jupyter notebook from rf pipeline: pipeline/eul1db/to_windows/_h/main.ipynb
__author__ = 'Michael Cuoco'

import pandas as pd
import sys,gc,traceback
from src.genome import Interval, Genome
    
def main():
    '''
    Convert eul1db_SRIPs.txt to a bed file, convert chr names to hs37d5 if necessary
    '''
    
    # load the eul1db file
    df0 = pd.read_csv(snakemake.input[0], sep="\t", skiprows=5, dtype={'chromosome':str})
    df0['chromosome'] = 'chr' + df0['chromosome']

    # filter for germline insertions from trusted studies
    df = df0[(df0['lineage']=='germline') & (df0['study_id'].isin(['Ewing2010', 'Ewing2011', 'Stewart2011', 'Beck2010', 'Brouha2002', 'Iskow2010'])) & (df0['g_start']==df0['g_stop'])]

    bed = df[['chromosome', 'g_start', 'g_stop']]
    bed = bed.rename(columns={'chromosome': 'chr', 'g_start': 'start', 'g_stop': 'end'})
    bed['start'] -= 1

    # fix names if necessary
    if snakemake.params.ref == 'hs37d5':
        # read in the chromosome map
        chrom_map = pd.read_csv("resources/hs37d5_map.tsv", sep="\t", names=["hs37d5", "ann"])

        for name in chrom_map["ann"].to_list():
            bed.loc[bed["chr"] == name, "chr"] = chrom_map.loc[ chrom_map["ann"] == name, "hs37d5"].values[0]

    # Filter out random chrs (3 instances)?
    # bed[~bed['chr'].str.contains("random")]
    bed.to_csv(snakemake.output[0], sep="\t", header=False, index=False)

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

    sys.stdout.close()