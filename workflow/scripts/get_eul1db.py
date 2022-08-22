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

    eul1db_pos = set()

    for (_, chrom, start, end) in df[['chromosome', 'g_start', 'g_stop']].itertuples():
        eul1db_pos.update([Interval(chrom, start, end)])
    
        len(eul1db_pos)

    genome = Genome(snakemake.input[0])
    xx = list(ig.windows_overlapping_intervals(genome, eul1db_pos, 750, 250))
    eul1 = pd.DataFrame.from_records((x.as_tuple() for x in xx), columns=['chrom', 'start', 'end']).set_index(['chrom', 'start', 'end'])
    eul1['in_NRdb'] = True # NR = non-reference
    eul1.to_csv(snakemake.output[0]) # why not use bed file for this?

if __name__ == '__main__':
    main()