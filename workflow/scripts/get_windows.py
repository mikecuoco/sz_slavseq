#!/usr/bin/env python
# adapted from Ricardo's jupyter notebook from rf pipeline: pipeline/eul1db/to_windows/_h/main.ipynb
__author__ = 'Michael Cuoco, Rohini Gadde'

import pandas as pd
from pyslavseq.genome import Interval, Genome
import pyslavseq.genome.interval_generator as ig
import sys, gc, traceback

def main():
    df = pd.read_csv(snakemake.input.srip[0], sep="\t", names=["chr", "start", "end"])

    db_pos = set()

    for (_, chrom, start, end) in df[['chr', 'start', 'end']].itertuples():
        db_pos.update([Interval(chrom, start, end)])
    
        len(db_pos)

    genome = Genome(snakemake.input.genome[0])
    xx = list(ig.windows_overlapping_intervals(genome, db_pos, 750, 250))
    l1 = pd.DataFrame.from_records((x.as_tuple() for x in xx), columns=['chr', 'start', 'end']).set_index(['chr', 'start', 'end'])
    l1['in_NRdb'] = True # NR = non-reference
    l1.to_csv(snakemake.output[0])

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