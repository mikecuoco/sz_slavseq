#!/usr/bin/env python
__author__ = "Apuã Paquola"


from .Interval import Interval, intersection


def _read_chromsizes_dict(chromsizes_fn):
    """Reads in a dict with chrom and size from a tab-delimited file"""
    f = open(chromsizes_fn)
    d = dict()
    for line in f:
        a = line.rstrip("\r\n").split("\t")
        if not a[1].isdigit():
            continue
        d[a[0]] = int(a[1])
    f.close()
    return d


class Genome:
    """A collection of intervals forming a genome"""

    def __init__(self, chromsizes_fn):
        self.chromsizes = _read_chromsizes_dict(chromsizes_fn)

    def fit_interval(self, iv):
        """If the interval iv is out of genome bounds, returns the intersection of it with the genome bounds"""
        if iv.chrom not in self.chromsizes:
            return None
        else:
            ci = Interval(iv.chrom, 0, self.chromsizes[iv.chrom])
            return intersection(iv, ci)
