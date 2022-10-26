#!/usr/bin/env python
__author__ = "Apuã Paquola"

import numpy as np
from . import Interval, overlaps


def chromosome_intervals(genome):
    """Yields one interval per each chromosome, spanning the whole chromosome.
    """
    for k in sorted(genome.chromsizes.keys()):
        yield Interval(k, 0, genome.chromsizes[k])


def windows_in_interval(genome, interval, window_size, step):
    """Yields regularly-spaced windows inside an interval.
    """
    for pos in range(interval.start, interval.end, step):
        yield genome.fit_interval(Interval(interval.chrom,
                                           pos,
                                           pos + window_size))


def windows_in_interval_collection(genome, collection, window_size, step):
    """Yields regularly-spaced windows inside intervals in a collection.
    """
    for interval in collection:
        for iv in windows_in_interval(genome, interval, window_size, step):
            yield iv


def windows_in_genome(genome, window_size, step):
    for iv in windows_in_interval_collection(genome, 
                                             chromosome_intervals(genome),
                                             window_size,
                                             step):
        yield iv


def hg19_dummy_intervals():
    for x in [('chr22', 25232000, 25432000), ('chrX', 129425250, 129625250), ('chr8', 72299000, 72499000), ('chr7', 62392750, 62592750), ('chr5', 157767250, 157967250), ('chr14', 60475250, 60675250), ('chr17', 46405000, 46605000), ('chr2', 46096000, 46296000), ('chr2', 42924000, 43124000), ('chr5', 35464750, 35664750), ('chr21', 1126250, 1326250), ('chr16', 24447000, 24647000), ('chr4', 23507000, 23707000)]:
        yield Interval(x[0], x[1], x[2])
    

def windows_centered_on_interval_centers(genome, intervals, window_size, step):
    a = set()
    for iv in intervals:
        newstart = round((iv.start + iv.end - window_size) / 2 / step) * step
        a.add(genome.fit_interval(
            Interval(iv.chrom, newstart, newstart+window_size)))
   
    for iv in sorted(a):
        yield iv


def windows_overlapping_intervals(genome, intervals, window_size, step):
    a = set()
    for iv in intervals:
        ni = Interval(iv.chrom,
                      ((iv.start - window_size) // step) * step,
                      ((iv.end + window_size + step) // step) * step)
        for wi in windows_in_interval(genome, ni, window_size, step):
            if wi is not None:
                if overlaps(iv, wi):
                    a.add(wi)
        
    for iv in sorted(a):
        yield iv


def random_empty_interval(genome):
    """ Returns a random interval in the genome, with start=end
    """
    k = sorted(genome.chromsizes.keys())
    p0 = np.array([genome.chromsizes[i] for i in k])

    chrom = np.random.choice(k, p=p0 / np.sum(p0))
    start = np.random.choice(genome.chromsizes[chrom])
    return Interval(chrom, start, start)
