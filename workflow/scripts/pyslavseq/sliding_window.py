#!/usr/bin/env python
# Created on: Jul 1, 2023 at 9:50:48 AM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)  # configure logging
from collections import deque
from itertools import chain, tee
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
    and r.mapping_quality >= 20
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
        minreads: int = 1,
        minrpm: float = 0,
    ) -> None:
        # save inputs to attributes
        self.bam = bam
        assert minreads > 0, "minreads must be > 0"

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
        rpmreads = ceil(minrpm * self.size_factor)
        self.minreads = max(minreads, rpmreads)

        # get all reads from bam
        self.reads = self.bam.fetch()
        self.reads = filter(self.read_filter, self.reads)
        self.reads = map(read_to_namedtuple, self.reads)

    def windows(
        self,
        reads: Generator,
        size: int = 200,
        step: int = 1,
    ) -> Generator:
        "Slide window across contig and yield each window with >= minreads reads"

        try:
            r = next(reads)
        except StopIteration:
            return

        w = deque()  # initialize window for reads
        n_dedup = 0  # initialize number of deduplicated reads
        count = 0  # initialize count of windows
        reflen = self.reflen  # store reference length
        minreads = self.minreads  # store minimum reads
        starts = range(0, reflen + 1 - size, step)
        ends = range(size, reflen + 1, step)
        for start, end in zip(starts, ends):
            # while w is not empty and the first read is outside the window
            while w and (w[0].three_end < start):
                n_dedup -= not w.popleft().is_duplicate

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
                            "reads": w,
                            "n_dedup": n_dedup,
                            "width": end - start,
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
                    "reads": w.copy(),
                    "n_dedup": n_dedup,
                    "width": end - start,
                }

    def findbg(
        self, windows: Generator, bg_windows: Generator, bgsize: int = 1e4
    ) -> Generator:
        """
        Find bg_windows centered around windows
        :param windows: generator of windows
        :param bg_windows: generator of background windows
        :param bgsize: size of background window
        """

        half_bgsize = bgsize / 2
        reflen_half_bgsize = self.reflen - half_bgsize

        for w in windows:
            center = int(w["Start"] + (w["width"] / 2))
            # find the background window that contains the center of the window
            # import pdb; pdb.set_trace()
            for bgw in bg_windows:
                # if the center of the window is at the beginning of the chromosome
                if center < half_bgsize:
                    break
                elif center == (bgw["Start"] + half_bgsize):
                    break
                # if the center of the window is at the end of the chromosome
                elif center > reflen_half_bgsize:
                    continue

            yield w, bgw

    def testbg(self, windows: Generator, bg_windows: Generator, mfold: int = 5):
        """
        Test if window signal is enriched above local background
        :param windows: generator of windows
        :param bg_windows: generator of background windows
        :param mfold: minimum fold change to consider window enriched above local background, default = 5
        """
        count = 0
        for w, bgw in self.findbg(windows, bg_windows):
            # if window is enriched above local background, yield it
            if ((w["n_dedup"] * bgw["width"]) / (bgw["n_dedup"] * w["width"])) > mfold:
                count += 1
                yield w

    def merge(
        self, windows: Generator, bandwidth: int = -100, refine: bool = False
    ) -> Generator:
        """
        Merge overlapping windows.
        :param windows: generator of windows
        :param bandwidth: maximum distance between windows to merge
        """

        windows = map(
            lambda x: {"Start": x["Start"], "End": x["End"], "reads": set(x["reads"])},
            windows,
        )

        try:
            p = next(windows)  # grab first window
        except StopIteration:
            return

        count = 0
        for n in windows:
            # if the next window is within the bandwidth, merge
            if n["Start"] <= (p["End"] + bandwidth):
                p["End"] = n["End"]
                p["reads"].update(n["reads"])
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
                p["width"] = p["End"] - p["Start"]
                yield p
                # start new window
                p = n

        count += 1
        p["width"] = p["End"] - p["Start"]
        yield p

        logger.info(f"Merged windows into {count} peaks on {self.contig}")

    def greedy_windows(
        self, windows: Generator, blacklist_size: int = 1000
    ) -> Generator:
        """
        Greedily choose windows with maximum number of reads
        :param windows: generator of windows
        :param blacklist_size: size of region to exclude around chosen window
        """

        # TODO: what to do with ties?
        # exhaust windows generator and sort by number of deduplicated reads
        windows = sorted(windows, key=lambda w: w["n_dedup"], reverse=True)
        blacklist = []  # initialize list of regions to exclude
        for w in windows:
            # if the window is not in the blacklist, yield it
            if not any(
                ((w["Start"] <= x[1]) and (w["End"] >= x[0])) for x in blacklist
            ):
                yield w
                blacklist.append(
                    (w["Start"] - blacklist_size, w["End"] + blacklist_size)
                )

        logger.info(
            f"Found {len(blacklist)} {self.size} bp windows using greedy algorithm"
        )

    def make_regions(
        self,
        mode: str = "peaks",
        bgtest: bool = False,
        bgsize: int = int(1e4),
        mfold: int = 5,
        greedy: bool = False,
        greedy_blacklist_size: int = 1000,
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
            logger.info(f"PROCESSING {c}")
            self.contig = c
            self.reflen = self.bam.get_reference_length(c)

            # get reads on contig
            chr_reads, self.reads = before_and_after(
                lambda r: r.reference_name == c, self.reads
            )

            # copy iterator for initial windows, background test, and background feature collection
            read_dict = {
                k: reads
                for k, reads in zip(
                    ["windows", "bgtest", 5e3, 1e4, 2e4], tee(chr_reads, 5)
                )
            }

            # define window generator
            windows = self.windows(reads=read_dict["windows"], **kwargs)

            if bgtest:
                bgkwargs = kwargs.copy()
                bgkwargs["size"] = int(bgsize)
                bg_windows = self.windows(reads=read_dict["bgtest"], **bgkwargs)
                windows = self.testbg(windows, bg_windows, mfold=mfold)

            if greedy:
                windows = self.greedy_windows(
                    windows, blacklist_size=greedy_blacklist_size
                )
                windows = sorted(
                    windows, key=lambda w: w["Start"]
                )  # re-sort windows by start position

            if mode == "peaks":
                windows = self.merge(windows, refine=refine)

            # add flanking windows for background features
            bgw = {}
            for bgsize, w in zip([5e3, 1e4, 2e4], tee(windows, 3)):
                bgkwargs = kwargs.copy()
                bgkwargs["size"] = int(bgsize)
                bg_windows = self.windows(reads=read_dict[bgsize], **bgkwargs)
                bgw[bgsize] = self.findbg(w, bg_windows, bgsize=bgsize)

            # yield windows
            for (w, bgw5), (_, bgw10), (_, bgw20) in zip(bgw[5e3], bgw[1e4], bgw[2e4]):
                w["Chromosome"] = c
                if collect_features:
                    w = features(w)
                    w["size_factor"] = self.size_factor
                    for bgsize, bg in zip([5e3, 1e4, 2e4], [bgw5, bgw10, bgw20]):
                        bg["Chromosome"] = c
                        for k, v in features(bg).items():
                            w[f"bg{int(bgsize)}_{k}"] = v
                    yield w
                else:
                    w["n_reads"] = len(w["reads"])
                    del w["reads"]
                    yield w
