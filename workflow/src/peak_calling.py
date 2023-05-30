#!/usr/bin/env python
# Created on: May 27, 2023 at 6:38:01 PM
__author__ = "Michael Cuoco"

# My custom peak callers
import pysam, re
import numpy as np
from scipy.stats import poisson
from collections import deque, namedtuple
from itertools import groupby, tee

# TODO: write tests

# TODO: classify peaks as ref or non-ref based on read2 alignment

# custom tuple to store reads in deque
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
        "mate_mapping_quality",
        "mate_cigar_string",
        "mate_reference_name",
        "mate_reference_start",
        "is_proper_pair",
    ],
)


class BasePeakCaller(object):
    def __init__(self, bam: pysam.AlignmentFile) -> None:
        self.bam = bam

    def read_converter(self, read: pysam.AlignedSegment) -> Read:
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
            read.get_tag("MQ") if read.has_tag("MQ") else None,
            read.get_tag("MC") if read.has_tag("MC") else None,
            read.next_reference_name,
            read.next_reference_start,
            read.is_proper_pair,
        )

    def read_filter(self, min_mapq=0) -> bool:
        "return callable to filter reads"
        return (
            lambda x: x.is_mapped
            and (not (x.is_secondary or x.is_supplementary))
            and (x.mapping_quality >= min_mapq)
        )


class SlidingPeakCaller(BasePeakCaller):
    def __init__(
        self,
        bam: pysam.AlignmentFile,
    ) -> None:
        super().__init__(bam)

    def window(
        self, reads, contig, bandwidth: int, step_size: int = 1, min_reads: int = 1
    ):
        """
        Generator that slides window of size bandwidth across a contig and yields windows with > 0 read alignments.
        :params reads: generator of reads
        :param contig: str, contig in bam
        :param window_size: int, size of window
        :param step_size: int, size of step bewteen windows
        :param min_reads: int, minimum number of reads in window
        """
        assert contig in self.bam.references, "contig {} not in bam".format(contig)

        # subfunction for window returning
        def window_dict():
            start = w[0].reference_start
            return {
                "Chromosome": r.reference_name,
                "Start": start,
                "End": r.reference_end,
                "width": r.reference_end - start,
                "center": (r.reference_end + start) // 2,
                "nreads": len(w),
            }

        # grab the first read
        try:
            r = next(reads)
        except StopIteration:
            return

        # initialize deque for sliding window
        # data structure for fast addition/removal of reads to/from sliding window
        w = deque()

        # iterate across the contig
        reflen = self.bam.get_reference_length(contig)
        for start in range(0, reflen - bandwidth + 1, step_size):
            end = start + bandwidth if start + bandwidth < reflen else reflen

            # remove reads outside window
            while w and w[0].reference_start < start:
                w.popleft()

            while (r.reference_start < end) and (r.reference_start >= start):
                w.append(r)
                try:
                    r = next(reads)
                except StopIteration:
                    if len(w):
                        yield window_dict()
                    else:
                        return

            # yield the window
            if len(w) > min_reads:
                yield window_dict()

    def test(self, w, lamb, pval_cutoff=1e-5):
        """
        Generator that yields a window from read alignments to contig using poisson test against background.
        :param w: window of reads from bam file
        :param lamb: float, background rate of reads (calculate from effective genome size)
        """

        # calculate p-value
        w["p"] = poisson._pmf(w["nreads"], lamb)
        w["fold_change"] = (
            np.float16(w["nreads"] / lamb) if lamb > 0 else np.float16(w["nreads"])
        )

        # yield window
        if w["p"] < pval_cutoff:
            return w

    def merge(self, windows):
        """
        Merge overlapping windows.
        :param windows: generator of windows
        """
        # merge the peaks
        try:
            pw = next(windows)  # grab first peak
        except StopIteration:  # skip if no peaks
            return

        for nw in windows:
            if nw["Start"] <= pw["End"]:
                pw["End"] = nw["End"]
                pw["width"] = pw["End"] - pw["Start"]
                pw["center"] = (pw["End"] + pw["Start"]) // 2
                pw["nreads"] += nw["nreads"]
            else:
                yield pw
                pw = nw

    def run(self, frag_size, bandwidth, eff_gen_len):
        "run the peak caller"

        reads = [r for r in self.get_reads()]
        nreads = len(reads)
        bg_lambda = (nreads * frag_size) / eff_gen_len

        # split reads by group and contig
        rg1, rg2 = tee(reads)
        f1 = lambda x: (x.is_read1 and x.is_forward) or (x.is_read2 and x.is_reverse)
        f2 = lambda x: (x.is_read1 and x.is_reverse) or (x.is_read2 and x.is_forward)

        for rg, f in zip([rg1, rg2], [f1, f2]):
            rg = filter(f, rg)
            rg = sorted(rg, lambda x: x.reference_name)
            rg = groupby(rg, lambda x: x.reference_name)
            for contig, reads in rg:
                windows = [
                    self.test(w, bg_lambda)
                    for w in self.window(reads, contig, bandwidth=bandwidth)
                ]
                merged = self.merge(windows)


class OverlapPeakCaller(BasePeakCaller):
    def __init__(self, bam: pysam.AlignmentFile) -> None:
        """
        Initializes a PeakCaller object with a bam file and peak and background window sizes.
        input:
        :param bam: pysam.AlignmentFile, opened bam file
        """
        super().__init__(bam)

    def cluster(self, reads, bandwidth: int = 0) -> deque:
        """
        Iterate over reads, cluster overlapping reads together
        """

        # define functions to classify reads as ref or nonref
        def parse_cigar(cigar):
            "convert cigar string to list of tuples of (length, operation)"
            return re.findall(r"(\d+)([MIDNSHP=X])", cigar)

        def ref_or_nonref(read, ref_reads, nonref_reads):
            "add to ref or nonref count depending on # clipped bases at start of read"
            cigar = parse_cigar(read.mate_cigar_string)

            # get cigar at start of read, accounting for strand
            end = cigar[-1] if read.is_forward else cigar[0]
            clipped = int(end[0]) if end[1] == "H" or end[1] == "S" else 0

            if read.is_proper_pair and clipped < 30:
                ref_reads += 1
            else:
                nonref_reads += 1
            return ref_reads, nonref_reads

        last_read = None
        ref_reads, nonref_reads = 0, 0
        peak = deque()

        for read in reads:
            if last_read is None:
                last_read = read
                if read.is_read1:
                    ref_reads, nonref_reads = ref_or_nonref(
                        read, ref_reads, nonref_reads
                    )
                peak.append(read)
                continue

            # if read is not overlapping with last read, yield as peak
            if read.reference_start > (last_read.reference_end + bandwidth):
                yield {
                    "Chromosome": read.reference_name,
                    "Start": peak[0].reference_start,
                    "End": last_read.reference_end,
                    "Strand": "-"
                    if (read.is_read1 and read.is_forward)
                    or (read.is_read2 and read.is_reverse)
                    else "+",
                    "nreads": len(peak),
                    "ref_reads": ref_reads,
                    "nonref_reads": nonref_reads,
                }
                # start new peak
                peak = deque()
                ref_reads, nonref_reads = 0, 0

            last_read = read
            ref_reads, nonref_reads = ref_or_nonref(read, ref_reads, nonref_reads)
            peak.append(read)

    def merge_peaks(self, peaks, bandwidth=0) -> dict:
        """
        Merge peaks that are within bandwidth of each other for a single contig
        """

        # initialize peak groups dict
        peak_groups = {
            "read1_fwd": None,
            "read2_fwd": None,
            "read1_rev": None,
            "read2_rev": None,
        }

        # iterate over peaks, merging paired peaks and yielding paired and unpaired
        for p in peaks:
            rg = p["read_group"]

            # if peak is unpaired, yield and overwrite
            if peak_groups[rg] is not None:
                yield {
                    "Chromosome": peak_groups[rg]["Chromosome"],
                    "Start": peak_groups[rg]["Start"],
                    "End": peak_groups[rg]["End"],
                    "Strand": "-" if rg == "read1_fwd" or "read2_rev" else "+",
                    "nreads": peak_groups[rg]["nreads"],
                }
            peak_groups[rg] = p

            # try to find paired peak
            for rg1, rg2 in zip(
                ["read1_fwd", "read1_rev", "read2_rev", "read2_fwd"],
                ["read2_rev", "read2_fwd", "read1_fwd", "read1_rev"],
            ):
                if rg == rg1 and (peak_groups[rg2] is not None):
                    # if concordant peak is within bandwidth of last peak, merge and yield
                    if peak_groups[rg1]["Start"] < (
                        peak_groups[rg2]["End"] + bandwidth
                    ):
                        yield {
                            "Chromosome": peak_groups[rg1]["Chromosome"],
                            "Start": peak_groups[rg2]["Start"],
                            "End": peak_groups[rg1]["End"],
                            "Strand": "-" if rg1 == "read1_fwd" or "read2_rev" else "+",
                            "nreads": peak_groups[rg1]["nreads"]
                            + peak_groups[rg2]["nreads"],
                        }
                        # reset
                        peak_groups[rg1] = None
                        peak_groups[rg2] = None

                    # return unmatched peaks
                    else:
                        yield {
                            "Chromosome": peak_groups[rg2]["Chromosome"],
                            "Start": peak_groups[rg2]["Start"],
                            "End": peak_groups[rg2]["End"],
                            "Strand": "-" if rg2 == "read1_fwd" or "read2_rev" else "+",
                            "nreads": peak_groups[rg2]["nreads"],
                        }
                        # reset
                        peak_groups[rg2] = None

        for p in peaks:
            rg = p["read_group"]

            # if peak is unpaired, yield and overwrite
            if peak_groups[rg] is not None:
                yield {
                    "Chromosome": peak_groups[rg]["Chromosome"],
                    "Start": peak_groups[rg]["Start"],
                    "End": peak_groups[rg]["End"],
                    "Strand": "-" if rg == "read1_fwd" or "read2_rev" else "+",
                    "nreads": peak_groups[rg]["nreads"],
                }
            peak_groups[rg] = p

            # try to find paired peak
            for rg1, rg2 in zip(
                ["read1_fwd", "read1_rev", "read2_rev", "read2_fwd"],
                ["read2_rev", "read2_fwd", "read1_fwd", "read1_rev"],
            ):
                if rg == rg1 and (peak_groups[rg2] is not None):
                    # if concordant peak is within bandwidth of last peak, merge and yield
                    if peak_groups[rg1]["Start"] < (
                        peak_groups[rg2]["End"] + bandwidth
                    ):
                        yield {
                            "Chromosome": peak_groups[rg1]["Chromosome"],
                            "Start": peak_groups[rg2]["Start"],
                            "End": peak_groups[rg1]["End"],
                            "Strand": "-" if rg1 == "read1_fwd" or "read2_rev" else "+",
                            "nreads": peak_groups[rg1]["nreads"]
                            + peak_groups[rg2]["nreads"],
                        }
                        # reset
                        peak_groups[rg1] = None
                        peak_groups[rg2] = None

                    # return unmatched peaks
                    else:
                        yield {
                            "Chromosome": peak_groups[rg2]["Chromosome"],
                            "Start": peak_groups[rg2]["Start"],
                            "End": peak_groups[rg2]["End"],
                            "Strand": "-" if rg2 == "read1_fwd" or "read2_rev" else "+",
                            "nreads": peak_groups[rg2]["nreads"],
                        }
                        # reset
                        peak_groups[rg2] = None

    def run(self, contig=None, bandwidth=0, min_mapq=0):
        """
        Run the peak caller.
        1. Get reads from bam
        2. Make peaks
        3. Merge peaks
        """

        # setup read group filters
        f1 = lambda x: (x.is_read1 and x.is_forward) or (x.is_read2 and x.is_reverse)
        f2 = lambda x: (x.is_read1 and x.is_reverse) or (x.is_read2 and x.is_forward)

        peaks = []
        merged = []

        for contig in self.bam.references:
            for f in [f1, f2]:
                # get the reads
                reads = filter(self.read_filter(min_mapq), self.bam.fetch(contig))
                reads = filter(f, reads)
                reads = map(self.read_converter, reads)

                # find peaks
                peaks += [p for p in self.cluster(reads, bandwidth=bandwidth)]

                # classify peaks

        # # get the peaks
        # peaks = [p for p in self.make_peaks(reads, r1_bandwidth)]

        # peaks = sorted(peaks, key=lambda x: (x["Chromosome"], x["Start"], x["End"]))
        # # sort peaks, split by contig
        # peaks = {c: list(g) for c, g in groupby(peaks, lambda x: x["Chromosome"])}

        # # merge the peaks
        # merged = [
        #     p for c in peaks.keys() for p in self.merge_peaks(peaks[c], r1_bandwidth)
        # ]
        # merged = sorted(merged, key=lambda x: (x["Chromosome"], x["Start"], x["End"]))

        # return peaks, merged
