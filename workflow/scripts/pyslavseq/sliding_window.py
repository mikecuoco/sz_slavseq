#!/usr/bin/env python
# Created on: Jul 1, 2023 at 9:50:48 AM
__author__ = "Michael Cuoco"

import re, logging, time
from collections import deque
from typing import Generator
from pysam import AlignmentFile, AlignedSegment
import numpy as np
import pandas as pd
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt
from .schemas import TAGS, FEATURES_SCHEMA, Read
import pyarrow as pa
import pyarrow.parquet as pq


AUTOSOMES = [f"chr{c}" for c in range(1, 23)]


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


class SlidingWindow(object):
    def __init__(
        self, bam: AlignmentFile, contigs: list = AUTOSOMES, min_mapq: int = 0
    ) -> None:
        self.bam = bam
        self.contigs = contigs
        self.min_mapq = min_mapq
        self.read_filter = (
            lambda x: x.is_read1
            and x.is_mapped
            and (not x.is_secondary)
            and (not x.is_supplementary)
            and (not x.is_duplicate)
        )

        total_reads = bam.count(read_callback=self.read_filter)
        self.size_factor = total_reads / 1e6
        logging.info(f"{total_reads} filtered reads in the bam file")

        # reset read filter to include min_mapq
        self.read_filter = (
            lambda x: x.is_read1
            and x.is_mapped
            and (not x.is_secondary)
            and (not x.is_supplementary)
            and (not x.is_duplicate)
            and (x.mapping_quality >= min_mapq)
        )

    def windows(self, reads, size: int = 200, step: int = 1) -> Generator:
        "Slide window across contig, yield windows with > min_rpm"
        try:
            r = next(reads)
            assert type(r) == Read, "Reads must be of type Read"
            contig = r.reference_name
            reflen = self.bam.get_reference_length(r.reference_name)
        except StopIteration:
            return

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

            yield {
                "Chromosome": w[0].reference_name,
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

        def peak_stats(p: dict):
            "calculate peak stats"

            p["n_reads"] = len(p["reads"])
            p["rpm"] = p["n_reads"] / self.size_factor

            # iterate over reads to calculate remaining stats
            p["n_ref_reads"], p["max_mapq"] = 0, 0
            starts = []
            for r in p["reads"]:
                p["max_mapq"] = max(p["max_mapq"], r.mapping_quality)
                p["n_ref_reads"] += r.isref_read
                if r.is_forward:
                    starts.append(r.reference_start)
                else:
                    starts.append(r.reference_end)
            p["n_unique_starts"] = len(set(starts))

            return p

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
                # calculate some stats
                p = peak_stats(p)
                # drop reads
                yield p
                # start new window
                assert (
                    p["Chromosome"] == n["Chromosome"]
                ), "Windows are not sorted by contig"
                p = n

    def features(self, w: dict) -> dict:
        """
        Extract features from a window of reads"
        :param w: window of reads
        """

        # initialize lists for features
        l = dict()
        for k in TAGS + ["3end", "5end", "mapq"]:
            l[k] = []
        f = FEATURES_SCHEMA.copy()  # initialize feature dictionary

        # add coordinates
        f["Chromosome"] = w["Chromosome"]
        f["Start"] = w["Start"]
        f["End"] = w["End"]

        # collect features from the reads in the window
        for r in w["reads"]:
            assert r.is_read1, "Reads must be all read 1"

            if r.is_forward:
                l["3end"].append(r.reference_end)
                l["5end"].append(r.reference_start)
            else:
                l["3end"].append(r.reference_start)
                l["5end"].append(r.reference_end)

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

        f["3end_gini"] = gini(np.array(l["3end"], dtype=np.float64))
        f["5end_gini"] = gini(np.array(l["5end"], dtype=np.float64))
        f["max_mapq"] = max(l["mapq"])

        # compute mean and quantiles for these features
        for tag in TAGS:
            if len(l[tag]) > 0:
                quantiles = np.quantile(l[tag], [0, 0.25, 0.5, 0.75, 1])
                for n, q in zip([0, 0.25, 0.5, 0.75, 1], quantiles):
                    f[tag + "_q" + str(n)] = q
                f[tag + "_mean"] = np.mean(l[tag])

        f["n_reads"] = f["n_fwd"] + f["n_rev"]
        f["rpm"] = f["n_reads"] / self.size_factor
        f["orientation_bias"] = np.maximum(f["n_fwd"], f["n_rev"]) / f["n_reads"]
        f["frac_proper_pairs"] = f["n_proper_pairs"] / f["n_reads"]

        # ensure output matches schema
        assert set(f.keys()) == set(FEATURES_SCHEMA.keys())

        return f

    def make_windows(
        self,
        strand_split: bool = False,
        merge: bool = False,
        features: bool = False,
        **kwargs,
    ) -> Generator:
        """
        Make windows on the given contigs
        :param strand_split: optionally generate windows separately for each strand
        :param merge: optionally merge overlapping windows
        :param features: optionally extract features from windows
        """

        # define read group filters
        if strand_split:
            rg = [
                lambda x: x.is_read1 and x.is_forward,
                lambda x: x.is_read1 and x.is_reverse,
            ]
        else:
            rg = [lambda x: x.is_read1]

        for c in self.contigs:
            logging.info(f"Making windows on {c}")
            for f in rg:
                reads = filter(self.read_filter, self.bam.fetch(c))
                reads = filter(f, reads)
                reads = map(read_to_namedtuple, reads)

                # define window generator
                if merge:
                    windows = self.merge(self.windows(reads, **kwargs))
                else:
                    windows = self.windows(reads, **kwargs)

                # yield windows
                for w in windows:
                    if features:
                        yield self.features(w)
                    else:
                        if "reads" in w:
                            del w["reads"]
                        yield w

    def write_windows(
        self, outfile: str, schema: pa.Schema, batch_size: int = 1000, **kwargs
    ) -> None:
        """
        Generate windows and write to disk
        :param outfile: path to output file
        :param schema: dictionary of column names and types, to be used as parquet schema
        """

        schema = pa.Schema.from_pandas(pd.Series(schema).to_frame().T)

        with pq.ParquetWriter(outfile, schema, compression="gzip") as writer:
            windows = []
            start = time.perf_counter()

            # write windows to disk in batches
            # TODO: make batches chromosomes
            for i, w in enumerate(self.make_windows(**kwargs)):
                windows.append(w)
                if (i > 0) and (
                    i % batch_size == 0
                ):  # write to disk every 1000 windows
                    writer.write_table(
                        pa.Table.from_pandas(pd.DataFrame(windows), schema=schema)
                    )
                    logging.info(
                        f"Processed {batch_size} windows in {time.perf_counter() - start:.2f} seconds"
                    )
                    windows = []
                    start = time.perf_counter()

            # write remaining windows
            if len(windows) > 0:
                writer.write_table(
                    pa.Table.from_pandas(pd.DataFrame(windows), schema=schema)
                )
                logging.info(
                    f"Processed {len(windows)} windows in {time.perf_counter() - start:.2f} seconds"
                )

    def coverage(
        self, peaks: pd.DataFrame, line1: pd.DataFrame, blacklist=None
    ) -> tuple:
        """
        Calculate coverage of peaks and annotations
        :param peaks: dataframe of peaks
        :param line1: dataframe of line1 annotations
        :param blacklist: dataframe of blacklist regions
        """

        # validate inputs
        for name, df in zip(["peaks", "knrgl", "rmsk"], [peaks, line1]):
            assert type(df) == pd.DataFrame, f"{name} must be a pandas dataframe"
            for c in ["Chromosome", "Start", "End"]:
                assert c in df.columns, f"{c} must be in {name} columns"
        assert "Name" in line1.columns, "Name must be in line1 columns"

        # ony use chromosomes of interest
        peaks = pr.PyRanges(peaks[peaks["Chromosome"].isin(self.contigs)])  # type: ignore
        line1 = pr.PyRanges(line1[line1["Chromosome"].isin(self.contigs)])  # type: ignore

        # optionally remove blacklist regions
        if blacklist is not None:
            blacklist = pr.PyRanges(
                blacklist[blacklist["Chromosome"].isin(self.contigs)]
            )
            peaks = peaks.overlap(blacklist, invert=True)  # type: ignore
            line1 = line1.overlap(blacklist, invert=True)  # type: ignore

        # label line1 by peaks
        line_df = line1.join(peaks, how="left").df

        # label peaks by line1
        peak_df = peaks.join(line1, how="left").sort().df  # type: ignore
        peak_df.loc[
            (peak_df.Name == "-1") & (peak_df.n_ref_reads == 0), "Name"
        ] = "NoneNR"
        peak_df.loc[
            (peak_df.Name == "-1") & (peak_df.n_ref_reads > 0), "Name"
        ] = "NoneR"
        peak_df["diff"] = peak_df.Start.diff()
        peak_df["width"] = peak_df.End - peak_df.Start

        # define order of L1s
        order = ["KNRGL", "L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]

        # PLOTTING
        # make subplots
        fig, axes = plt.subplots(2, 3, figsize=(13, 8))
        fig.tight_layout(pad=3)

        # plot 1: how many peak reads overlap each L1?
        sns.ecdfplot(
            line_df,
            x="n_reads",
            hue="Name",
            hue_order=order,
            ax=axes[0, 0],
        ).set(xscale="log", xlim=(1, None), title="Peak reads in L1s")

        # plot 2: how many reads are in each peak?
        sns.ecdfplot(
            peak_df,
            x="n_reads",
            hue="Name",
            hue_order=order + ["NoneR", "NoneNR"],
            ax=axes[0, 1],
        ).set(xscale="log", xlim=(1, None), title="Reads in Peaks by label")

        # plot 3: how wide are the peaks?
        sns.ecdfplot(
            peak_df,
            x="width",
            hue="Name",
            hue_order=order + ["NoneR", "NoneNR"],
            ax=axes[0, 2],
        ).set(title="Peak Width by label")

        # plot 4: how are nrefreads distributed?
        sns.ecdfplot(
            peak_df,
            x="n_ref_reads",
            hue="Name",
            hue_order=order + ["NoneR", "NoneNR"],
            ax=axes[1, 0],
        ).set(xlim=(1, None), xscale="log", title="Reference reads in Peaks by label")

        peak_df["frac_ref_reads"] = peak_df.n_ref_reads / peak_df.n_reads
        sns.ecdfplot(
            peak_df,
            x="frac_ref_reads",
            hue="Name",
            hue_order=order + ["NoneR", "NoneNR"],
            ax=axes[1, 1],
        ).set(title="Fraction of reference reads in Peaks by label")

        # plot 6: how many peaks labelled by each class?
        axes[1, 2] = peak_df.groupby(["Name"]).size().plot.barh(rot=0)
        axes[1, 2].bar_label(axes[1, 2].containers[0])
        axes[1, 2].set_xscale("log")
        axes[1, 2].set_title("Peak Count by label")
        axes[1, 2].set_xlim(1, max(axes[1, 2].get_xlim()) * 2)

        return axes, peak_df, line_df
