#!/usr/bin/env python
# Created on: Jul 1, 2023 at 9:50:48 AM
__author__ = "Michael Cuoco"

import logging, time

logger = logging.getLogger(__name__)  # configure logging
from collections import deque
from typing import Generator
from pysam import AlignmentFile
import numpy as np
from math import ceil
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from .features import features, basic_stats, read_to_namedtuple

CHROMOSOMES = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]


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
                w.popleft()

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
        r = next(gen_regions)  # get first region
        regions = []  # initialize list of regions

        schema = pa.Schema.from_pandas(pd.DataFrame([r]))
        schema = schema.append(pa.field("rpm", pa.float64()))

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
                    logger.info(f"Skipping {c} because it is empty")
                    continue

                # convert regions to dataframe and remove empty regions
                regions = pd.DataFrame(regions).query("n_reads > 0")

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
                try:
                    writer.write_table(pa.Table.from_pandas(regions, schema=schema))
                except:
                    import pdb

                    pdb.set_trace()
                # reset regions list
                regions = []
