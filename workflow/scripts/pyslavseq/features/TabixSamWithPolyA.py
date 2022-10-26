#!/usr/bin/env python
__author__ = "Apuã Paquola"

from .TabixSam import TabixSam


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


class TabixSamWithPolyA(TabixSam):
    """ Extracts length of maximal poly A sequence in read2 as a feature from a TabixSam
    """

    def _r2_poly_N_length(self, s, letter, weight):
        ws = [weight if s[i] == letter else -(1 - weight) for i in range(len(s))]
        (maxsum, start, end) = maxsumseq(ws)
        return end - start

    def fetch(self, chrom, start, end):
        for r in TabixSam.fetch(self, chrom, start, end):
            if 'Y2' in r['tags']:
                r['r2_poly_a_length'] = self._r2_poly_N_length(r['tags']['Y2'], 'A', 0.2)
            yield r
