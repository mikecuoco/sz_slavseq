#!/usr/bin/env python
__author__ = "Apuã Paquola"

import argparse
import pysam
from ..genome import interval_generator as ig
from ..genome import Genome, Interval
from .WindowFeatures import WindowFeatures
from .TabixSamWithPolyA import TabixSamWithPolyA


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--min_mapq', type=float, required=True)
    parser.add_argument('--min_ya', type=float, required=True)
    parser.add_argument('--max_yg', type=float, required=True)
    parser.add_argument('--chromsizes', required=True)
    parser.add_argument('--window_size', type=int, required=True)
    parser.add_argument('--window_step', type=int, required=True)
    parser.add_argument('--min_secondary_mapq', type=float, required=True)
    parser.add_argument('--write_header_to')
    parser.add_argument('--chrom')
    parser.add_argument('--start', type=int)
    parser.add_argument('--end', type=int)
    parser.add_argument('--pilot', default=False, action='store_true')
    # parser.add_argument('--alupoly', default=False, action='store_true')
    parser.add_argument('filename')
    args = parser.parse_args()

    tabixsam = TabixSamWithPolyA(pysam.Tabixfile(args.filename))

    genome = Genome(args.chromsizes)

    if args.pilot:
        wg_iter = ig.windows_in_interval_collection(genome, ig.hg19_dummy_intervals(), args.window_size,
                                                    args.window_step)
    elif args.chrom is not None:
        wg_iter = ig.windows_in_interval(genome, Interval(args.chrom, args.start, args.end), args.window_size,
                                         args.window_step)
    else:
        wg_iter = ig.windows_in_genome(genome, args.window_size, args.window_step)

    k = None
    # i=0
    for w in wg_iter:
        # i+=1
        # if i>1000:
        #    break

        wf = WindowFeatures(tabixsam, w.chrom, w.start, w.end, args.min_mapq, args.min_ya, args.max_yg,
                            args.min_secondary_mapq)
        f = wf.features()
        if k is None:
            k = sorted(f.keys())

            if args.write_header_to is not None:
                of = open(args.write_header_to, 'w')
            else:
                of = None
            print('\t'.join(['chrom', 'start', 'end'] + k), file=of)

            if args.write_header_to is not None:
                of.close()

        print('\t'.join([str(x) for x in [w.chrom, w.start, w.end]] + [str(f[x]) for x in k]))


if __name__ == "__main__":
    main()
