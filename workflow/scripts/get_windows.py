#!/usr/bin/env python

import pandas as pd
from pyslavseq.genome import Interval, Genome
import pyslavseq.genome.interval_generator as ig

def main():
    db_pos = set()

    for (_, chrom, start, end) in df[['chr', 'start', 'stop']].itertuples():
        db_pos.update([Interval(chrom, start, end)])
    
        len(db_pos)

    genome = Genome(snakemake.input[0])
    xx = list(ig.windows_overlapping_intervals(genome, db_pos, 750, 250))
    l1 = pd.DataFrame.from_records((x.as_tuple() for x in xx), columns=['chr', 'start', 'end']).set_index(['chr', 'start', 'end'])
    l1['in_NRdb'] = True # NR = non-reference
    l1.to_csv(snakemake.output[0])

if __name__ == '__main__':
    main()