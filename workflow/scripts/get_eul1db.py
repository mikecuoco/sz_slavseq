#!/usr/bin/env python
# adapted from Ricardo's jupyter notebook from rf pipeline: pipeline/eul1db/to_windows/_h/main.ipynb
__author__ = 'Michael Cuoco'

import pandas as pd
import sys,gc,traceback
from pyslavseq.genome import Interval, Genome
import pyslavseq.genome.interval_generator as ig
    
def main():
    # load the eul1db file
    df0 = pd.read_csv(snakemake.input.eul1db, sep="\t", skiprows=5, dtype={'chromosome':str})
    df0['chromosome'] = 'chr' + df0['chromosome']

    # filter for germline insertions from trusted studies
    df = df0[(df0['lineage']=='germline') & (df0['study_id'].isin(['Ewing2010', 'Ewing2011', 'Stewart2011', 'Beck2010', 'Brouha2002', 'Iskow2010'])) & (df0['g_start']==df0['g_stop'])]

    # generate intervals for db
    eul1db_pos = set()
    for (_, chrom, start, end) in df[['chromosome', 'g_start', 'g_stop']].itertuples():
        eul1db_pos.update([Interval(chrom, start, end)])

    # generate intervals for genome
    genome = Genome(snakemake.input.genome[0])

    # place eul1db intervals into genome
    xx = list(ig.windows_overlapping_intervals(genome, eul1db_pos, 750, 250))
    eul1 = pd.DataFrame.from_records((x.as_tuple() for x in xx), columns=['chrom', 'start', 'end']).set_index(['chrom', 'start', 'end'])
    eul1['in_NRdb'] = True # NR = non-reference
    eul1.to_csv(snakemake.output[0]) # why not use bed file for this?

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