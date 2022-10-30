#!/usr/bin/env python
__author__ = "Apuã Paquola"

import re

def maxsumseq(sequence):
    start, end, sum_start = -1, -1, -1
    maxsum_, sum_ = 0, 0
    for i, x in enumerate(sequence):
        sum_ += x
        if maxsum_ < sum_:  # found maximal subsequence so far
            maxsum_ = sum_
            start, end = sum_start, i
        elif sum_ < 0:  # start new sequence
            sum_ = 0
            sum_start = i
    # assert maxsum_ == maxsum(sequence)
    assert maxsum_ == sum(sequence[start + 1:end + 1])
    return maxsum_, start + 1, end + 1


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
        
    def _r2_poly_N_length(self, s, letter, weight):
        ws = [weight if s[i] == letter else -(1 - weight) for i in range(len(s))]
        (maxsum, start, end) = maxsumseq(ws)
        return end - start

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

    def fetch_from_file(self, chrom, start, end):
        return self._set_primary_r1(self._fetch(chrom, start, end))

    # Extracts the length of maximal poly A sequence in read2 as a feature from a TabixSam
    def fetch_polyA(self, chrom, start, end):
        for r in self.fetch_from_file(chrom, start, end):
            if 'Y2' in r['tags']:
                r['r2_poly_a_length'] = self._r2_poly_N_length(r['tags']['Y2'], 'A', 0.2)
            yield r
