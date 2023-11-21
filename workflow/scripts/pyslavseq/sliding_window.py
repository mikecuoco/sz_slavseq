#!/usr/bin/env python
# Created on: Jul 1, 2023 at 9:50:48 AM
__author__ = "Michael Cuoco"

import re, logging, time

logger = logging.getLogger(__name__)  # configure logging
from collections import deque
from typing import Generator
from pysam import AlignmentFile, AlignedSegment
import numpy as np
from math import ceil
import pandas as pd
from .schemas import TAGS, Read
import pyarrow as pa
import pyarrow.parquet as pq

CHROMOSOMES = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]


def isref_read(read: AlignedSegment) -> bool:
    "return True if read is ref, False if non-ref based on mate tag"
    assert read.is_read1, "Read must be read 1"

    # get cigar at start of read, accounting for strand
    if not read.has_tag("MC"):
        return False

    cigar = re.findall(r"(\d+)([MIDNSHP=X])", str(read.get_tag("MC")))
    end = cigar[-1] if not read.is_reverse else cigar[0]
    clipped = int(end[0]) if end[1] == "H" or end[1] == "S" else 0

    return read.is_proper_pair and (clipped < 30)


def read_to_namedtuple(read: AlignedSegment) -> Read:
    "convert pysam.AlignedSegment to hashable namedtuple Read"
    return Read(
        read.reference_name,
        read.reference_start,
        read.reference_end,
        read.infer_read_length(),
        read.is_read1,
        read.is_read2,
        not read.is_reverse,
        read.is_reverse,
        read.mapping_quality,
        len(read.get_tag("SA").split(";")[:-1]) if read.has_tag("SA") else 0,  # type: ignore
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
        read.is_duplicate,
    )


def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq: http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    array = array.flatten()  # all values are treated equally, arrays must be 1d
    if np.amin(array) < 0:
        array -= np.amin(array)  # values cannot be negative
    array += 0.0000001  # values cannot be 0
    array = np.sort(array)  # values must be sorted
    index = np.arange(1, array.shape[0] + 1)  # index per array element
    n = array.shape[0]  # number of array elements
    return (np.sum((2 * index - n - 1) * array)) / (
        n * np.sum(array)
    )  # Gini coefficient


def basic_stats(p: dict):
    """
    Calculate basic stats about a region
    :param p: dict of region, must contain "reads" key
    """

    for n in ["Chromosome", "Start", "End", "reads"]:
        assert n in p, f"Region must contain {n}"

    # iterate over reads to calculate stats
    p["n_ref_reads"], p["max_mapq"], p["n_duplicates"] = 0, 0, 0
    starts = []
    fwd, rev = 0, 0
    for r in p["reads"]:
        if r.is_duplicate:
            p["n_duplicates"] += 1
            continue

        if r.is_reverse:
            rev += 1
        else:
            fwd += 1
        p["max_mapq"] = max(p["max_mapq"], r.mapping_quality)
        p["n_ref_reads"] += r.isref_read

        if r.is_forward:
            starts.append(r.reference_start)
        else:
            starts.append(r.reference_end)
    p["n_unique_starts"] = len(set(starts))
    p["n_reads"] = fwd + rev

    # add strand info
    if (fwd == 0) and (rev > 0):
        p["Strand"] = "-"
    elif (rev == 0) and (fwd > 0):
        p["Strand"] = "+"

    # remove reads
    del p["reads"]

    return p


def features(p: dict) -> dict:
    """
    Extract features from a window of reads
    All reads must be read1
    :param p: dict of region, must contain "reads" key
    """

    for n in ["Chromosome", "Start", "End", "reads"]:
        assert n in p, f"Region must contain {n}"

    # initialize lists for features
    l = dict()
    for k in TAGS + ["3end", "5end", "mapq", "starts"]:
        l[k] = []

    # initialize features dict
    f = {
        "Chromosome": p["Chromosome"],
        "Start": p["Start"],
        "End": p["End"],
        "n_reads": 0,
        "n_fwd": 0,
        "n_rev": 0,
        "n_duplicates": 0,
        "n_proper_pairs": 0,
        "n_ref_reads": 0,
        "3end_gini": float(0),
        "5end_gini": float(0),
        "max_mapq": 0,
        "rpm": float(0),
        "orientation_bias": float(0),
        "frac_proper_pairs": float(0),
        "n_unique_starts": 0,
    }

    for tag in TAGS:
        for n in [0, 0.25, 0.5, 0.75, 1]:
            f[tag + "_q" + str(n)] = float(0)
        f[tag + "_mean"] = float(0)

    # collect features from the reads in the window
    for i, r in enumerate(p["reads"]):
        if not r.is_read1:
            raise Exception("Reads must all be read1")

        if i == 0:
            start = r.reference_start

        if r.is_duplicate:
            f["n_duplicates"] += 1
            continue

        l["starts"].append(r.reference_start)
        if r.is_forward:
            l["3end"].append(r.reference_end - start)
            l["5end"].append(r.reference_start - start)
        else:
            l["3end"].append(r.reference_start - start)
            l["5end"].append(r.reference_end - start)

        l["mapq"].append(r.mapping_quality)
        f["n_proper_pairs"] += r.is_proper_pair
        f["n_ref_reads"] += r.isref_read

        if r.is_reverse:
            f["n_rev"] += 1
        else:
            f["n_fwd"] += 1

        for tag in TAGS:
            if "_normed" in tag:
                if getattr(r, tag.replace("_normed", "")):
                    if tag in [
                        "L1_alignment_score_normed",
                        "mate_alignment_score_normed",
                    ]:  # adjust alignments scores for read length
                        l[tag].append(
                            getattr(r, tag.replace("_normed", ""))
                            / getattr(r, "mate_read_length")
                        )
                    elif tag in ["alignment_score_normed"]:
                        l[tag].append(
                            getattr(r, tag.replace("_normed", ""))
                            / getattr(r, "read_length")
                        )
            elif getattr(r, tag):
                l[tag].append(getattr(r, tag))

    # compute mean and quantiles for these features
    for tag in TAGS:
        if len(l[tag]) > 0:
            quantiles = np.quantile(l[tag], [0, 0.25, 0.5, 0.75, 1])
            for n, q in zip([0, 0.25, 0.5, 0.75, 1], quantiles):
                f[tag + "_q" + str(n)] = q
            f[tag + "_mean"] = np.mean(l[tag])

    f["n_reads"] = f["n_fwd"] + f["n_rev"]

    # if reads are empty, return empty features
    if f["n_reads"] == 0:
        return f

    f["3end_gini"] = gini(np.array(l["3end"], dtype=np.float64))
    f["5end_gini"] = gini(np.array(l["5end"], dtype=np.float64))
    f["max_mapq"] = max(l["mapq"])
    f["n_unique_starts"] = len(set(l["starts"]))

    return f


class SlidingWindow(object):
    def __init__(
        self,
        bam: AlignmentFile,
        contigs: list | None = None,
        read_filter: callable = lambda x: x.is_read1,
        mode: str = "peaks",
        collect_features: bool = False,
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
        self.mode = mode
        self.collect_features = collect_features

        # count reads in bam satisfying read filter
        total_reads = bam.count(read_callback=self.read_filter)
        self.size_factor = total_reads / 1e6
        logger.info(f"{total_reads} filtered reads in the bam file")

    def windows(
        self,
        reads,
        size: int = 200,
        step: int = 1,
        min_reads: int = 0,
        min_rpm: float = 0,
    ) -> Generator:
        "Slide window across contig and yield each window"
        try:
            r = next(reads)
            assert type(r) == Read, "Reads must be of type Read"
            contig = r.reference_name
            reflen = self.bam.get_reference_length(r.reference_name)
        except StopIteration:
            return

        rpm_reads = ceil(min_rpm * self.size_factor)
        min_reads = max(min_reads, rpm_reads)
        w = deque()
        for start in range(0, reflen + 1, step):
            end = start + size if start + size < reflen else reflen

            # while w is not empty and the first read is outside the window
            while w and w[0].reference_start < start:
                w.popleft().isref_read

            # while the next read is inside the window
            while r.reference_start < end:
                w.append(r)
                try:
                    r = next(reads)
                    assert r.reference_name == contig, "Reads are not sorted by contig"
                except StopIteration:
                    break

            # return window if it has min_rpm
            if len(w) >= min_reads:
                yield {
                    "Chromosome": contig,
                    "Start": start,
                    "End": end,
                    "reads": set(w),
                }

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

        for n in windows:
            # if the next window is within the bandwidth, merge
            if n["Start"] <= (p["End"] + bandwidth):
                p["End"] = n["End"]
                p["reads"] = p["reads"].union(n["reads"])
            # otherwise, yield the previous window and start a new one
            else:
                # yield merged window
                yield p
                assert (
                    p["Chromosome"] == n["Chromosome"]
                ), "Windows are not sorted by contig"
                # start new window
                p = n

    def make_regions(self, **kwargs) -> Generator:
        """
        Make windows on the given contigs
        :param kwargs: arguments to pass to windows()
        """

        for c in self.contigs:
            reads = filter(self.read_filter, self.bam.fetch(c))
            reads = map(read_to_namedtuple, reads)

            # define window generator
            if self.mode == "peaks":
                windows = self.merge(self.windows(reads, **kwargs))
            elif self.mode == "windows":
                windows = self.windows(reads, **kwargs)
            else:
                raise Exception(
                    f"Mode {self.mode} is not valid, must be 'peaks' or 'windows'"
                )

            # yield windows
            for w in windows:
                if self.collect_features:
                    yield features(w)
                else:
                    yield basic_stats(w)

    def write_regions(self, outfile: str, **kwargs) -> None:
        """
        Generate regions and write to disk
        :param outfile: path to output file
        :param kwargs: arguments to pass to make_regions(), to pass to windows()
        """

        # create generator of regions and get first region
        gen_regions = self.make_regions(**kwargs)
        r = next(gen_regions)

        schema = pa.Schema.from_pandas(pd.Series(r).to_frame().T)

        regions = []  # initialize list of regions
        # write regions to disk in batches by chromosome
        with pq.ParquetWriter(outfile, schema, compression="gzip") as writer:
            for c in self.contigs:
                start = time.perf_counter()

                # collect regions for this chromosome
                while r["Chromosome"] == c:
                    regions.append(r)
                    try:
                        r = next(gen_regions)
                    except StopIteration:
                        break

                logger.info(
                    f"Generated {len(regions)} {self.mode} on {c} in {time.perf_counter() - start:.2f} seconds"
                )

                # skip empty chromosomes
                if len(regions) == 0:
                    continue

                # convert regions to dataframe and remove empty regions
                regions = pd.DataFrame(regions).query("n_reads > 0").fillna(0)

                # compute metrics on vector scale
                regions["rpm"] = regions["n_reads"] / self.size_factor
                if self.collect_features:
                    regions["orientation_bias"] = (
                        np.maximum(regions["n_fwd"], regions["n_rev"])
                        / regions["n_reads"]
                    )
                    regions["frac_proper_pairs"] = (
                        regions["n_proper_pairs"] / regions["n_reads"]
                    )
                    regions["frac_duplicates"] = regions["n_duplicates"] / (
                        regions["n_reads"] + regions["n_duplicates"]
                    )

                # write to disk
                logger.info(f"Writing {len(regions)} {self.mode} on {c} to disk")
                writer.write_table(pa.Table.from_pandas(regions, schema=schema))

                # reset regions list
                regions = []
