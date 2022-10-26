#!/usr/bin/env python
__author__ = "Apuã Paquola"

import itertools
import numpy as np
import re

# Utility functions


def cigar2reflen(cigar):
    r = 0
    for m in re.finditer('([0-9]+)[MDN]', cigar):
        r += int(m.group(1))
    return r


def start_pos(read):
    p = read['sam'][3]
    if read['sam'][1] & 16:
        p = p + cigar2reflen(read['sam'][5]) - 1
    return p


def end_pos(read):
    p = read['sam'][3]
    if (read['sam'][1] & 16) == 0:
        p = p + cigar2reflen(read['sam'][5]) - 1
    return p


def distance0(chr1, pos1, chr2, pos2):
    if chr1 != chr2:
        return np.inf
    return abs(pos2 - pos1)


def distance(read1, read2):
    if (read1['sam'][1] & 4) != 0 or (read2['sam'][1] & 4) != 0:
        return np.inf
    return distance0(read1['sam'][2], start_pos(read1), read2['sam'][2], start_pos(read2))


# Read list functions

def primary_read_filter(_iter, func=lambda x: True):
    """
    Only pass this filter the reads whose primary r1 aligment satisfy func().
    """

    l = list(_iter)

    selected_qnames = set()

    for x in l:
        if (x['primary_r1'] is x) and func(x):
            selected_qnames.add(x['sam'][0])

    for x in l:
        if x['sam'][0] in selected_qnames:
            yield x


def unmapped_r2_count(alignments):
    """ Return the # of read pairs whose r2 is unmapped """
    return len([1 for x in alignments if (x['sam'][1] & (128 + 4)) == (128 + 4)])


def secondary_read_count(alignments, min_secondary_mapq):
    """ Returns the number of primary r1 reads that have secondary alignments, either in r1 or r2.
    """
    return len(set(x['sam'][0] for x in alignments if (x['sam'][1] & 2304) != 0 and x['sam'][4] >= min_secondary_mapq))


def secondary_read_distances(alignments, min_secondary_mapq):
    """ Returns the number of primary r1 reads that have secondary alignments, either in r1 or r2.
    """
    return [distance(x['primary_r1'], x) for x in alignments if
            (x['sam'][1] & 2304) != 0 and x['sam'][4] >= min_secondary_mapq]


def r1r2_distances(alignments):
    return [
        distance(x['primary_r1'], x)
        for x in alignments
        if ((x['sam'][1] & 2304) == 0) and ((x['sam'][1] & 128) != 0)
    ]


def pretty_print_dict(d):
    for x in d:
        print(x, d[x], sep='\t')


class WindowFeatures:
    """ Given a tabixsam file, a window in the genome, and filtering parameters,
    returns a set of features for that window.
    """

    def __init__(self, tabixsam, chrom, start, end, min_mapq, min_ya, max_yg, min_secondary_mapq, library_3_or_5,
                 ensearch=None):
        self.tabixsam = tabixsam
        self.chrom = chrom
        self.start = start
        self.end = end
        self.min_mapq = min_mapq
        self.min_ya = min_ya
        self.max_yg = max_yg
        self.min_secondary_mapq = min_secondary_mapq
        self.library_3_or_5 = library_3_or_5
        self.ensearch = ensearch

    def r2_starts_and_strands_(self, alignments):
        # returns r2 strand
        # ignore r2 outside window
        return [
            (start_pos(x), [-1, 1][(x['sam'][1] & 16) == 0])
            for x in alignments
            if ((x['sam'][1] & 2304) == 0) and ((x['sam'][1] & 128) != 0)
               and x['sam'][2] == self.chrom and self.start <= x['sam'][3] <= self.end]

    def r1_ends_and_reverse_strands_(self, alignments):
        # returns r2 strand
        return [
            (end_pos(x), [-1, 1][(x['sam'][1] & 16) != 0])
            for x in alignments
            if ((x['sam'][1] & 2304) == 0) and ((x['sam'][1] & 64) != 0)
        ]

    def read_start_peak_features(self, positions_and_strands):
        if len(positions_and_strands) == 0:
            return np.nan, np.nan, np.nan

        read_starts = [x[0] for x in positions_and_strands]
        ifreq = np.array(np.unique(read_starts, return_counts=True)).transpose()
        imax = np.argmax(ifreq[:, 1])
        peakpos = ifreq[imax, 0]
        peakcount = ifreq[imax, 1]

        peakstrands = [x[1] for x in positions_and_strands if x[0] == peakpos]
        ifreq = np.array(np.unique(peakstrands, return_counts=True)).transpose()
        imax = np.argmax(ifreq[:, 1])
        peak_te_strand = ifreq[imax, 0]

        if self.library_3_or_5 == '5':
            peak_te_strand = -peak_te_strand

        return peakpos, peakcount, peak_te_strand

    def features0(self, alignments):
        primary_r1_reads = list(itertools.filterfalse(lambda x: not (x['primary_r1'] is x), alignments))

        feature_dict = dict()
        feature_dict['count'] = len(primary_r1_reads)
        feature_dict['plus_strand.count'] = len(
            list(itertools.filterfalse(lambda x: not (x['sam'][1] & 16 == 0), primary_r1_reads)))
        feature_dict['secondary.count'] = secondary_read_count(alignments, self.min_secondary_mapq)
        feature_dict['median_r2_poly_a_length'] = np.nan if feature_dict['count'] == 0 else np.median(
            [x['r2_poly_a_length'] for x in primary_r1_reads])

        feature_dict['median_ys'] = np.nan if feature_dict['count'] == 0 else np.median(
            [x['tags']['YS'] for x in primary_r1_reads])
        feature_dict['median_yg'] = np.nan if feature_dict['count'] == 0 else np.median(
            [x['tags']['YG'] for x in primary_r1_reads])
        feature_dict['median_ya'] = np.nan if feature_dict['count'] == 0 else np.median(
            [x['tags']['YA'] for x in primary_r1_reads])
        feature_dict['median_xd'] = np.nan if feature_dict['count'] == 0 else np.median(
            [x['tags']['XD'] for x in primary_r1_reads])

        d = r1r2_distances(alignments)
        feature_dict['median_r1r2_distance'] = np.nan if len(d) == 0 else np.median(d)
        feature_dict['r1r2_distance_1k_200k.count'] = len([x for x in d if 1000 <= x <= 200000])

        sd = secondary_read_distances(alignments, self.min_secondary_mapq)
        feature_dict['secondary.median_distance'] = np.nan if len(sd) == 0 else np.median(sd)

        feature_dict['secondary_distance_1k_200k.count'] = len([x for x in d if 1000 <= x <= 200000])

        feature_dict['unmapped_r2_count'] = unmapped_r2_count(alignments)

        positions_and_strands = self.r2_starts_and_strands_(alignments)
        (peakpos, peakcount, peak_te_strand) = self.read_start_peak_features(positions_and_strands)

        feature_dict['r2_peak_position'] = peakpos
        feature_dict['r2_peak_count'] = peakcount
        feature_dict['r2_peak_te_strand'] = peak_te_strand

        if self.ensearch is not None:
            (motif, en_pos, en_score) = self.ensearch.pos_and_score(self.chrom, peakpos, peak_te_strand)
            feature_dict['r2_en_motif'] = motif
            feature_dict['r2_en_position'] = en_pos
            feature_dict['r2_en_score'] = en_score

        positions_and_strands = self.r1_ends_and_reverse_strands_(alignments)
        (peakpos, peakcount, peak_te_strand) = self.read_start_peak_features(positions_and_strands)

        feature_dict['r1_peak_position'] = peakpos
        feature_dict['r1_peak_count'] = peakcount
        feature_dict['r1_peak_te_strand'] = peak_te_strand

        if self.ensearch is not None:
            (motif, en_pos, en_score) = self.ensearch.pos_and_score(self.chrom, peakpos, peak_te_strand)
            feature_dict['r1_en_motif'] = motif
            feature_dict['r1_en_position'] = en_pos
            feature_dict['r1_en_score'] = en_score

        return feature_dict

    def features(self):
        all_alignments = list(self.tabixsam.fetch(self.chrom, self.start, self.end))
        alignments_from_well_mapped_reads = list(
            primary_read_filter(all_alignments, lambda x: x['sam'][4] >= self.min_mapq))

        feature_dict = dict()

        f0 = self.features0(alignments_from_well_mapped_reads)
        for k in f0.keys():
            feature_dict['all_reads.' + k] = f0[k]

        alignments_from_well_mapped_yayg_reads = list(primary_read_filter(alignments_from_well_mapped_reads,
                                                                          lambda x: x['tags']['YA'] >= self.min_ya and
                                                                                    x['tags']['YG'] <= self.max_yg))

        f0 = self.features0(alignments_from_well_mapped_yayg_reads)
        for k in f0.keys():
            feature_dict['yayg_reads.' + k] = f0[k]

        alignments_from_well_mapped_yg_reads = list(
            primary_read_filter(alignments_from_well_mapped_reads, lambda x: x['tags']['YG'] <= self.max_yg))

        f0 = self.features0(alignments_from_well_mapped_yg_reads)
        for k in f0.keys():
            feature_dict['yg_reads.' + k] = f0[k]

        return feature_dict
