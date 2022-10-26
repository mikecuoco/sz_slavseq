#!/usr/bin/env python
__author__ = "Apuã Paquola"

import pysam
import MOODS.tools
import MOODS.scan

from ..genome import Interval
from ..util import rc
import numpy as np
import functools


class ENSearch:
    """ Search for a L1 endonuclease motif given by a PWM from Jurka 1997.
    """

    def __init__(self, genome, genome_fasta_file, left_flank, right_flank):
        # Matrix from PMID 9050872 Jurka 1997
        self.matrix = [[60, 71, 279, 248, 238, 241],
                       [34, 37, 3, 4, 9, 14],
                       [43, 26, 32, 72, 72, 46],
                       [207, 210, 30, 20, 25, 43]]

        self.genome = genome
        self.fa = pysam.Fastafile(genome_fasta_file)
        self.left_flank = left_flank
        self.right_flank = right_flank
        self.threshold = MOODS.tools.threshold_from_p(self.matrix, MOODS.tools.flat_bg(4), 0.2)

    @functools.lru_cache(maxsize=1024, typed=False)
    def pos_and_score(self, chrom, pos, te_strand):
        if pos is np.nan:
            return (np.nan, np.nan, np.nan)

        if te_strand == 1:
            # if te is in the + strand
            start = pos - self.left_flank
            end = pos + self.right_flank
        else:
            start = pos - self.right_flank
            end = pos + self.left_flank

        iv = self.genome.fit_interval(Interval(chrom, start, end))

        s = self.fa.fetch(iv.chrom, iv.start, iv.end).upper()

        if te_strand == 1:
            zeropos = pos - iv.start
        else:
            s = rc(s)
            zeropos = iv.end - pos
            # print("\n>> ", chrom, pos, te_strand,  file=sys.stderr, flush=True)

        # print(">> ", s,  file=sys.stderr, flush=True)

        def moods_results(s, matrix, left_flank):

            # This is a hack to work around a bug in MOODS. The last 3
            # parameters are: int iterations = 10, unsigned int MULT = 2,
            # size_t LIMIT_MULT = 10
            # Setting MULT and LIMIT_MULT to the size of the string fixes it

            results = MOODS.scan.scan_best_hits_dna(s, [matrix], 1, 10, len(s), len(s))
            if len(results) == 0:
                en_pos = - self.left_flank
                en_score = 0
                motif = ''
            else:
                # print(">> ", len(results[0]),  file=sys.stderr, flush=True)
                spos = results[0][0].pos
                motif = s[spos:spos + len(matrix[0])]
                en_pos = spos - zeropos
                en_score = int(results[0][0].score)

                return motif, en_pos, en_score

        return moods_results(s, self.matrix, self.left_flank)
