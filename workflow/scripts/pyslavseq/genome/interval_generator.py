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
