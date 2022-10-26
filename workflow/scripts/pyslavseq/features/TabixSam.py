#!/usr/bin/env python
__author__ = "Apuã Paquola"

import re


class TabixSam:
    """TabixSam is a file format in which the 3 first columns are chrom,start,end, and the rest are the SAM columns
    """

    def __init__(self, tabixfile):
        """tabixfile - an already opened tabix file.
        Example: tabixsam = TabixSam(pysam.Tabixfile(filename))
        """
        self.tabixfile = tabixfile

    def _fetch(self, chrom, start, end):
        try:
            # This try-except block is a workaround. Tabix fetch raises a ValueError: invalid region `b'chr1:1-750'`
            # if chr1 is not present in the index

            for line in self.tabixfile.fetch(chrom, start, end):
                a = line.split('\t')
                tags = dict()
                for f in a[14:]:
                    m = re.search('^(..):(.):(.*)', f)
                    if m.group(2) == 'i':
                        tags[m.group(1)] = int(m.group(3))
                    else:
                        tags[m.group(1)] = m.group(3)

                sam = [int(a[i]) if i in [4, 6, 7, 10, 11] else a[i] for i in range(3, 14)]
                yield {'main_pos': a[0:3],
                       'sam': sam,
                       'tags': tags}

        except (KeyError, ValueError):
            pass

    def _set_primary_r1(self, _iter):
        primary_r1 = dict()
        l = list(_iter)
        for x in l:
            x['primary_r1'] = None
            if (x['sam'][1] & 64) and ((x['sam'][1] & 2304) == 0):
                primary_r1[x['sam'][0]] = x
        for x in l:
            x['primary_r1'] = primary_r1[x['sam'][0]]
            yield x

    def fetch(self, chrom, start, end):
        return self._set_primary_r1(self._fetch(chrom, start, end))
