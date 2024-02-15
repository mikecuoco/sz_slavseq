#!/usr/bin/env python
# Created on: Jul 1, 2023 at 9:50:48 AM
__author__ = "Michael Cuoco"

import logging, time

logger = logging.getLogger(__name__)  # configure logging
from collections import deque
from itertools import chain, tee, takewhile
from typing import Generator
from pysam import AlignmentFile
from math import ceil
import numpy as np
from .features import features, read_to_namedtuple

CHROMOSOMES = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]

read_filter = (
    lambda r: (not r.is_read2)
    and r.is_mapped
    and (not r.is_supplementary)
    and (not r.is_secondary)
)


def before_and_after(predicate, it):
    """Variant of takewhile() that allows complete
    access to the remainder of the iterator.

    >>> it = iter('ABCdEfGhI')
    >>> all_upper, remainder = before_and_after(str.isupper, it)
    >>> ''.join(all_upper)
    'ABC'
    >>> ''.join(remainder)     # takewhile() would lose the 'd'
    'dEfGhI'

    Note that the true iterator must be fully consumed
    before the remainder iterator can generate valid results.
    """
    transition = []

    def true_iterator():
        for elem in it:
            if predicate(elem):
                yield elem
            else:
                transition.append(elem)
                return

    return true_iterator(), chain(transition, it)


class SlidingWindow(object):
    def __init__(
        self,
        bam: AlignmentFile,
        contigs: list | None = None,
        read_filter: callable = read_filter,
    ) -> None:
        # save inputs to attributes
        self.bam = bam

        # validate contigs
        if contigs:
            for c in contigs:
                if c not in CHROMOSOMES:
                    raise Exception(f"{c} is not valid chromosome")
            self.contigs = contigs
        else:
            self.contigs = []
            for i in range(bam.nreferences):
                if bam.get_reference_name(i) in CHROMOSOMES:
                    self.contigs.append(bam.get_reference_name(i))

        self.read_filter = read_filter

        # get library size
        libsize_read_filter = (
            lambda r: (not r.is_read2)
            and r.is_mapped
            and (not r.is_supplementary)
            and (not r.is_secondary)
            and (not r.is_duplicate)
        )
        self.library_size = bam.count(read_callback=libsize_read_filter)
        self.size_factor = self.library_size / 1e6
        logger.info(f"{self.library_size} filtered reads in the bam file")

        # get all reads from bam
        self.reads = self.bam.fetch()
        self.reads = filter(self.read_filter, self.reads)
        self.reads = map(read_to_namedtuple, self.reads)

    def windows(
        self,
        reads: Generator,
        size: int = 200,
        step: int = 1,
        minreads: int = 0,
        minrpm: float = 0,
    ) -> Generator:
        "Slide window across contig and yield each window"

        try:
            r = next(reads)
        except StopIteration:
            return

        rpmreads = ceil(minrpm * self.size_factor)
        minreads = max(minreads, rpmreads)
        w = deque(maxlen=self.library_size)  # initialize window for reads
        n_dedup = 0  # initialize number of deduplicated reads
        count = 0  # initialize count of windows
        reflen = self.reflen  # store reference length
        for start in range(0, reflen + 1, step):
            end = start + size
            end = end if end < reflen else reflen

            # while w is not empty and the first read is outside the window
            while w and w[0].three_end < start:
                r_out = w.popleft()
                n_dedup -= not r_out.is_duplicate

            # while the next read is inside the window
            while r.three_end < end:
                w.append(r)
                n_dedup += not r.is_duplicate
                try:
                    r = next(reads)
                except StopIteration:
                    if n_dedup >= minreads:
                        count += 1
                        yield {
                            "Start": start,
                            "End": end,
                            "reads": set(w),
                        }
                    logger.info(
                        f"Generated {count} {size} bp windows on {self.contig} with >= {minreads} reads"
                    )
                    return

            # return window if it has minreads
            if n_dedup >= minreads:
                count += 1
                yield {
                    "Start": start,
                    "End": end,
                    "reads": set(w),
                }

    def testbg(
        self, windows: Generator, bg_windows: Generator, mfold: int = 5
    ) -> Generator:
        """
        Test if window signal is enriched above local background
        :param windows: generator of windows
        :param bg_windows: generator of background windows
        :param mfold: minimum fold change to consider window enriched above local background, default = 5
        """

        windows = filter(lambda x: len(x["reads"]) > 0, windows)
        w = next(windows)
        wsize = w["End"] - w["Start"]
        half_wsize = wsize / 2
        windows = chain([w], windows)

        bgw = next(bg_windows)
        bgsize = bgw["End"] - bgw["Start"]
        bg_windows = chain([bgw], bg_windows)

        half_bgsize = bgsize / 2
        reflen_half_bgsize = self.reflen - half_bgsize
        count = 0
        for w in windows:
            center = w["Start"] + half_wsize
            # find the background window that contains the center of the window
            for bgw in bg_windows:
                # if the center of the window is at the beginning of the chromosome
                if center < half_bgsize:
                    break
                if (bgw["Start"] + half_bgsize) == center:
                    break
                # if the center of the window is at the end of the chromosome
                if center > reflen_half_bgsize:
                    continue

            # test for fold enrichment
            if ((len(w["reads"]) / wsize) / (len(bgw["reads"]) / bgsize)) > mfold:
                count += 1
                yield w

        logger.info(
            f"Found {count} {wsize} bp windows {mfold}x enriched above {bgsize} bp background on {self.contig}"
        )

    def merge(
        self, windows: Generator, bandwidth: int = -100, refine: bool = False
    ) -> Generator:
        """
        Merge overlapping windows.
        :param windows: generator of windows
        :param bandwidth: maximum distance between windows to merge
        """

        # filter for non-empty windows
        windows = filter(lambda x: len(x["reads"]) > 0, windows)

        p = next(windows)  # grab first window

        count = 0
        for n in windows:
            # if the next window is within the bandwidth, merge
            if n["Start"] <= (p["End"] + bandwidth):
                p["End"] = n["End"]
                p["reads"] = p["reads"].union(n["reads"])
            # otherwise, yield the previous window and start a new one
            else:
                # yield merged window
                # refine start and end first
                if refine:
                    start, end = p["End"], p["Start"]
                    for r in p["reads"]:
                        if r.three_end < start:
                            start = r.three_end
                        if r.three_end > end:
                            end = r.three_end
                    p["Start"], p["End"] = start, end
                count += 1
                yield p
                # start new window
                p = n

        logger.info(f"Merged windows into {count} peaks on {self.contig}")

    def make_regions(
        self,
        mode: str = "peaks",
        bgsize: int = int(1e4),
        bgtest: bool = False,
        mfold: int = 5,
        refine: bool = False,
        collect_features: bool = False,
        **kwargs,
    ) -> Generator:
        """
        Make windows on the given contigs
        :param mode: peaks or windows. Default = "peaks"
        :param bgtest: whether to test for fold enrichment of windows above local background. Default = False
        :param bgsize: size of local background window to test against. Only applied if bgtest = True. Default = 10000
        :param collect_features: whether to collect features for the regions. Default = False
        :param kwargs: arguments to pass to windows()
        """
        if mode not in ["windows", "peaks"]:
            raise Exception(f"Mode {mode} is not valid, must be 'peaks' or 'windows'")

        for c in self.contigs:
            self.contig = c
            self.reflen = self.bam.get_reference_length(c)

            # get reads on contig
            reads, self.reads = before_and_after(
                lambda r: r.reference_name == c, self.reads
            )
            if bgtest:
                bgreads, reads = tee(reads)

            # define window generator
            windows = self.windows(reads=reads, **kwargs)

            # TODO: handle if testbg results in no windows
            if bgtest:
                bgkwargs = kwargs.copy()
                bgkwargs["size"] = bgsize
                bg_windows = self.windows(reads=bgreads, **bgkwargs)
                windows = self.testbg(windows, bg_windows, mfold=mfold)

            if mode == "peaks":
                windows = self.merge(windows, refine=refine)

            # yield windows
            for w in windows:
                w["Chromosome"] = c
                if collect_features:
                    w = features(w)
                    w["size_factor"] = self.size_factor
                    yield w
                else:
                    w["n_reads"] = len(w["reads"])
                    del w["reads"]
                    yield w
