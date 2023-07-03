#!/usr/bin/env python
# Created on: Jul 1, 2023 at 9:50:48 AM
__author__ = "Michael Cuoco"

import re, logging
from collections import deque, namedtuple
from typing import Generator
from pysam import AlignmentFile, AlignedSegment
import pandas as pd

Read = namedtuple(
    "Read",
    [
        "reference_name",
        "reference_start",
        "reference_end",
        "is_read1",
        "is_read2",
        "is_forward",
        "is_reverse",
        "mapping_quality",
        "cigar_string",
        "alignment_score",
        "L1_alignment_score",
        "L1_reference_start",
        "L1_reference_end",
        "L1_Acount",
        "mate_alignment_score",
        "mate_read_length",
        "mate_mapping_quality",
        "mate_cigar_string",
        "mate_reference_name",
        "mate_reference_start",
        "is_proper_pair",
        "isref_read",
    ],
)

AUTOSOMES = [f"chr{c}" for c in range(1, 22)]


def isref_read(read: AlignedSegment) -> bool:
    "return True if read is ref, False if non-ref based on mate tag"
    assert read.is_read1, "Read must be read 1"

    # get cigar at start of read, accounting for strand
    if not read.has_tag("MC"):
        return False

    cigar = re.findall(r"(\d+)([MIDNSHP=X])", read.get_tag("MC"))
    end = cigar[-1] if read.is_forward else cigar[0]
    clipped = int(end[0]) if end[1] == "H" or end[1] == "S" else 0

    return read.is_proper_pair and (clipped < 30)


def _read_to_namedtuple(read: AlignedSegment) -> Read:
    "convert pysam.AlignedSegment to namedtuple Read"
    return Read(
        read.reference_name,
        read.reference_start,
        read.reference_end,
        read.is_read1,
        read.is_read2,
        read.is_forward,
        read.is_reverse,
        read.mapping_quality,
        read.cigarstring,
        read.get_tag("AS") if read.has_tag("AS") else None,
        read.get_tag("L1") if read.has_tag("L1") else None,
        read.get_tag("LS") if read.has_tag("LS") else None,
        read.get_tag("LE") if read.has_tag("LE") else None,
        read.get_tag("LA") if read.has_tag("LA") else None,
        read.get_tag("MS") if read.has_tag("MS") else None,
        read.get_tag("ML") if read.has_tag("ML") else None,
        read.get_tag("MQ") if read.has_tag("MQ") else None,
        read.get_tag("MC") if read.has_tag("MC") else None,
        read.next_reference_name,
        read.next_reference_start,
        read.is_proper_pair,
        isref_read(read),
    )


class SlidingWindow(object):
    def __init__(
        self, bam: AlignmentFile, contigs: list = AUTOSOMES, min_mapq: int = 0
    ) -> None:
        self.bam = bam
        self.contigs = contigs
        self.min_mapq = min_mapq
        self.read_filter = (
            lambda x: x.is_mapped
            and (not (x.is_secondary or x.is_supplementary))
            and (x.mapping_quality >= min_mapq)
        )

    def windows(self, reads, size: int, step: int, min_reads: int = 1) -> Generator:
        "Slide window across contig, yield windows with >= min_reads"
        # TODO: change min_reads to min_RPM
        try:
            r = next(reads)
            assert type(r) == Read, "Reads must be of type Read"
            contig = r.reference_name
            reflen = self.bam.get_reference_length(r.reference_name)
        except StopIteration:
            return

        w = deque()
        for start in range(0, reflen - size + 1, step):
            end = start + size if start + size < reflen else reflen

            # while w is not empty and the first read is outside the window
            while w and w[0].reference_start < start:
                w.popleft()

            # while the next read is inside the window
            while r.reference_start < end:
                w.append(r)
                try:
                    r = next(reads)
                    assert r.reference_name == contig, "Reads are not sorted by contig"
                except StopIteration:
                    break

            # yield the window
            if len(w) >= min_reads:
                yield {
                    "Chromosome": w[0].reference_name,
                    "Start": w[0].reference_start,
                    "End": w[-1].reference_end,
                    "reads": set(w),
                }

    def merge(self, windows: Generator, bandwidth: int = 0) -> Generator:
        """
        Merge overlapping windows.
        :param windows: generator of windows
        :param bandwidth: maximum distance between windows to merge
        """

        try:
            p = next(windows)  # grab first window
        except StopIteration:  # skip if no peaks
            return

        for n in windows:
            # if the next window is within the bandwidth, merge
            if n["Start"] <= (p["End"] + bandwidth):
                p["End"] = n["End"]
                p["reads"] = p["reads"].union(n["reads"])
            # otherwise, yield the previous window and start a new one
            else:
                # calculate some window stats
                p["nreads"] = len(p["reads"])
                p["nrefreads"] = 0
                p["max_mapq"] = 0
                for r in p["reads"]:
                    p["max_mapq"] = max(p["max_mapq"], r.mapping_quality)
                    p["nrefreads"] += r.isref_read
                yield p
                # start new window
                assert (
                    p["Chromosome"] == n["Chromosome"]
                ), "Windows are not sorted by contig"
                p = n

    def call_peaks(self, **kwargs) -> pd.DataFrame:
        "Run the peak caller on the given contigs"

        # define read groups
        f1 = lambda x: (x.is_read1 and x.is_forward)
        f2 = lambda x: (x.is_read1 and x.is_reverse)

        peaks = []
        for c in self.contigs:
            logging.info(f"Calling peaks on {c}")
            for f in [f1, f2]:
                reads = filter(self.read_filter, self.bam.fetch(c))
                reads = filter(f, reads)
                reads = map(_read_to_namedtuple, reads)
                peaks += [r for r in self.merge(self.windows(reads, **kwargs))]

        return pd.DataFrame(peaks)

    def make_windows(self, **kwargs) -> Generator:
        "Make windows on the given contigs"

        for c in self.contigs:
            logging.info(f"Making windows on {c}")
            reads = filter(lambda x: x.is_read1, self.bam.fetch(c))
            reads = filter(self.read_filter, reads)
            reads = map(_read_to_namedtuple, reads)
            for w in self.windows(reads, **kwargs):
                yield w
