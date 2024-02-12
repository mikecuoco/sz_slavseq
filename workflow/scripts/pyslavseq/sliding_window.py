#!/usr/bin/env python
# Created on: Jul 1, 2023 at 9:50:48 AM
__author__ = "Michael Cuoco"

import logging, time

logger = logging.getLogger(__name__)  # configure logging
from collections import deque
from typing import Generator
from pysam import AlignmentFile
from math import ceil
import numpy as np
from scipy.stats import poisson, false_discovery_control
from .features import features, read_to_namedtuple

CHROMOSOMES = [f"chr{c}" for c in range(1, 23)]

read_filter = (
    lambda r: (not r.is_read2)
    and r.is_mapped
    and (not r.is_supplementary)
    and (not r.is_secondary)
)


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
        library_size = bam.count(read_callback=libsize_read_filter)
        self.size_factor = library_size / 1e6
        logger.info(f"{library_size} filtered reads in the bam file")

    def windows(
        self,
        size: int = 200,
        step: int = 1,
        minreads: int = 0,
        minrpm: float = 0,
    ) -> Generator:
        "Slide window across contig and yield each window"

        self.wsize = size

        # create generator of reads from contig from bam file
        reads = filter(
            self.read_filter, self.bam.fetch(self.contig, multiple_iterators=True)
        )
        reads = map(read_to_namedtuple, reads)

        try:
            r = next(reads)
        except StopIteration:
            return

        rpmreads = ceil(minrpm * self.size_factor)
        minreads = max(minreads, rpmreads)
        logger.info(
            f"Generating {size} bp windows on {self.contig} with minreads: {minreads}"
        )
        w = deque()  # initialize window for reads
        n_dedup = 0  # initialize number of deduplicated reads
        for start in range(0, self.reflen + 1, step):
            end = start + size if start + size < self.reflen else self.reflen

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
                        return {
                            "Start": start,
                            "End": end,
                            "reads": set(w),
                        }
                    else:
                        return

            # return window if it has minreads
            if n_dedup >= minreads:
                yield {
                    "Start": start,
                    "End": end,
                    "reads": set(w),
                }

    def bgtest(self, windows: Generator, bg_windows: Generator) -> Generator:
        """
        Test if window signal is significantly enriched above local background
        """

        # TODO: do fold enrichment test rather than p-value test
        cov, bgcov = np.array([]), np.array([])
        windows = list(windows)  # get values from generator into memory

        half_bgsize = self.bgsize / 2
        reflen_half_bgsize = self.reflen - (self.bgsize / 2)
        logger.info(
            f"Testing windows agaist {self.bgsize} local background on {self.contig}"
        )
        for w in windows:
            center = (w["Start"] + w["End"]) / 2
            for bgw in bg_windows:
                if center < half_bgsize:
                    break
                bg_center = (bgw["Start"] + bgw["End"]) / 2
                if bg_center == center:
                    break
                if center > reflen_half_bgsize:
                    continue

            cov = np.append(cov, len(w["reads"]))
            bgcov = np.append(bgcov, len(bgw["reads"]))

        cov = cov / self.wsize
        bgcov = bgcov / self.bgsize
        p_values = 1 - poisson._cdf(cov, bgcov)
        if self.bgqval:
            cutoff = self.bgqval
            p_values = false_discovery_control(p_values)  # BH correction
        else:
            cutoff = self.bgpval

        for w, p in zip(windows, p_values):
            if p < cutoff:
                yield w

    def merge(self, windows: Generator, bandwidth: int = 0) -> Generator:
        """
        Merge overlapping windows.
        :param windows: generator of windows
        :param bandwidth: maximum distance between windows to merge
        """

        # filter for non-empty windows
        windows = filter(lambda x: len(x["reads"]) > 0, windows)

        try:
            p = next(windows)  # grab first window
        except StopIteration:  # skip if no windows left
            return

        logger.info(f"Merging overlapping windows on {self.contig}")
        for n in windows:
            # if the next window is within the bandwidth, merge
            if n["Start"] <= (p["End"] + bandwidth):
                p["End"] = n["End"]
                p["reads"] = p["reads"].union(n["reads"])
            # otherwise, yield the previous window and start a new one
            else:
                # yield merged window
                # refine start and end first
                start, end = p["End"], p["Start"]
                for r in p["reads"]:
                    if r.three_end < start:
                        start = r.three_end
                    if r.three_end > end:
                        end = r.three_end
                p["Start"], p["End"] = start, end
                yield p
                # start new window
                p = n

    def make_regions(
        self,
        mode: str = "peaks",
        size: int = 200,
        bgsize: int = int(1e4),
        bgtest: bool = False,
        bgpval: float | None = 0.05,
        bgqval: float | None = None,
        collect_features: bool = False,
        **kwargs,
    ) -> Generator:
        """
        Make windows on the given contigs
        :param mode: peaks or windows. Default = "peaks"
        :param bgtest: whether to statistiacllyl test for enrichment of windows above local background. Default = False
        :param bgsize: size of local background window to test against. Only applied if bgtest = True. Default = 10000
        :param collect_features: whether to collect features for the regions. Default = False
        :param kwargs: arguments to pass to windows()
        """
        if mode not in ["windows", "peaks"]:
            raise Exception(f"Mode {mode} is not valid, must be 'peaks' or 'windows'")

        for c in self.contigs:
            self.contig = c
            self.reflen = self.bam.get_reference_length(c)
            windows = self.windows(size, **kwargs)

            if bgtest:
                self.bgsize = bgsize
                self.bgpval = bgpval
                self.bgqval = bgqval
                if bgpval and bgqval:
                    raise Exception("bgpval and bgqval cannot both be specified")
                bgkwargs = kwargs.copy()
                bg_windows = self.windows(size=bgsize, **bgkwargs)
                windows = self.bgtest(windows, bg_windows)

            # define window generator
            if mode == "peaks":
                windows = self.merge(windows, bandwidth=int(size / 2) * -1)

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
