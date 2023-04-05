#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import pysam
import pandas as pd
import numpy as np
from scipy.stats import entropy
from collections import deque, namedtuple, defaultdict
import sys


def reads(bam, contig):
    """Yield reads from a bam file for a given contig."""
    for r in bam.fetch(contig, multiple_iterators=True):
        if not r.is_duplicate and not r.is_unmapped:
            yield r


# define named tuple to hold window information
Window = namedtuple("Window", ["chr", "start", "end", "reads"])


def windows(bam, contig, window_size=750, step_size=250):
    """
    Yield windows of a given size and step size from a given contig.
    :param bam: pysam.AlignmentFile
    :param contig: str
    :param window_size: int
    :param step_size: int
    """

    # get the first read
    # if there are no reads, return
    try:
        read_iter = reads(bam, contig)
        r = next(read_iter)
    except StopIteration:
        return

    # use double-ended queue to hold reads
    read_queue = deque()
    reflen = bam.get_reference_length(contig)
    for start in range(0, reflen - window_size + 1, step_size):
        end = start + window_size if start + window_size < reflen else reflen

        # while read_queue is not empty and the first read is outside the window
        while read_queue and read_queue[0].reference_start < start:
            read_queue.popleft()

        # while the next read is inside the window
        while r.reference_start < end:
            read_queue.append(r)
            try:
                r = next(read_iter)
            except StopIteration:
                break

        yield Window(contig, start, end, read_queue)


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


def window_features(window):
    """Compute features of a window."""

    l = defaultdict(list)
    f = defaultdict(float)

    for r in window.reads:
        if r.is_read1:
            l["r1_mapq"].append(r.mapping_quality)
            l["r1_starts"].append(float(r.reference_start))
            l["template_lengths"].append(r.template_length)
            if r.is_reverse:
                f["r1_rev"] += 1
            else:
                f["r1_fwd"] += 1

            for tag in ["YG", "YA", "YS"]:
                if r.has_tag(tag):
                    l[tag].append(r.get_tag(tag))
                    f[tag + "_reads"] += 1

        else:
            l["r2_mapq"].append(r.mapping_quality)
            l["r2_starts"].append(float(r.reference_start))
            if r.is_reverse:
                f["r2_rev"] += 1
            else:
                f["r2_fwd"] += 1

    f["r1_reads"] = f["r1_fwd"] + f["r1_rev"]
    f["r2_reads"] = f["r2_fwd"] + f["r2_rev"]
    f["total_reads"] = f["r1_reads"] + f["r2_reads"]

    if f["r1_reads"] > 0:
        f["frac_r1_fwd"] = f["r1_fwd"] / f["r1_reads"]
        f["r1_orientation_entropy"] = entropy(
            [f["frac_r1_fwd"], 1 - f["frac_r1_fwd"]], base=2
        )
        f["r1_mean_mapq"] = np.mean(l["r1_mapq"])
        f["r1_sd_mapq"] = np.std(l["r1_mapq"])
        f["r1_starts_gini"] = gini(np.array(l["r1_starts"]))
        f["mean_template_length"] = np.mean(l["template_lengths"])
        f["sd_template_length"] = np.std(l["template_lengths"])
        for tag in ["YG", "YA", "YS"]:
            if f[tag + "_reads"] > 0:
                f[tag + "_mean"] = np.mean(l[tag])
        f["YA_YG_ratio"] = f["YA_mean"] / (f["YG_mean"] + 0.01)
    if f["r2_reads"] > 0:
        f["frac_r2_fwd"] = f["r2_fwd"] / f["r2_reads"]
        f["r2_orientation_entropy"] = entropy(
            [f["frac_r2_fwd"], 1 - f["frac_r2_fwd"]], base=2
        )
        f["r2_mean_mapq"] = np.mean(l["r2_mapq"])
        f["r2_sd_mapq"] = np.std(l["r2_mapq"])
        f["r2_starts_gini"] = gini(np.array(l["r2_starts"]))

    f["chr"] = window.chr
    f["start"] = window.start
    f["end"] = window.end

    return f


def window_filter(window):
    """Only yield windows that pass filter."""

    # if window["total_reads"] >= 3 and window["ya_mean"] > window["yg_mean"]:
    # if window["r1_fwd"] >= 3 or window["r2_fwd"] >= 3:
    if window["total_reads"] >= 1:
        return True
    else:
        return False


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    # open bam file
    bam = pysam.AlignmentFile(snakemake.input["bam"], "rb")

    # make windows and get features
    df = []

    WINDOW_SIZE = snakemake.params["window_size"]
    STEP_SIZE = snakemake.params["window_step"]
    BG_WINDOW_SIZE = int(snakemake.params["window_size"] * 10)
    BG_STEP_SIZE = int(snakemake.params["window_step"] / 10)

    for contig in bam.references:

        reflen = bam.get_reference_length(contig)
        # make generator for background windows (10kb)
        bg_windows = windows(
            bam, contig, window_size=BG_WINDOW_SIZE, step_size=BG_STEP_SIZE
        )

        # get first window
        try:
            bg = next(bg_windows)
            bg_center = bg.start + (bg.end - bg.start) / 2
        except StopIteration:
            continue

        # iterate over windows
        for w in windows(bam, contig, window_size=WINDOW_SIZE, step_size=STEP_SIZE):
            f = window_features(w)
            if not window_filter(f):
                continue

            w_center = w.start + (w.end - w.start) / 2
            while True:
                if (
                    (w_center == bg_center)
                    or (w_center < bg_center and bg.start == 0)
                    or (
                        w_center > bg_center
                        and bg.start > reflen - BG_WINDOW_SIZE - BG_STEP_SIZE
                    )
                ):
                    ff = window_features(bg)
                    for k in ff.keys():
                        f[k + "_bg"] = ff[k]
                    break
                else:
                    try:
                        bg = next(bg_windows)
                        bg_center = bg.start + (bg.end - bg.start) / 2
                    except StopIteration:
                        break

            # if "chr_bg" not in f:
            #     import pdb
            #     pdb.set_trace()

            assert (
                "chr_bg" in f
            ), "failed to find background window for window at {}:{}-{}".format(
                w.chr, w.start, w.end
            )

            df.append(f)

    df = pd.DataFrame(df).fillna(0)
    df.rename(
        {"chr": "Chromosome", "start": "Start", "end": "End"}, axis=1, inplace=True
    )
    df.drop(["chr_bg", "start_bg", "end_bg"], axis=1, inplace=True)

    # add cell_id and donor_id columns
    df["cell_id"] = snakemake.wildcards.sample
    df["donor_id"] = snakemake.wildcards.donor

    # save
    df.to_parquet(snakemake.output[0])

    sys.stderr.close()
