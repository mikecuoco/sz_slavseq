#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import pysam
import pandas as pd
from collections import deque, namedtuple
import sys


def reads(bam, contig):
    """Yield reads from a bam file for a given contig."""
    for r in bam.fetch(contig, multiple_iterators=True):
        if not r.is_duplicate:
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


def window_features(window):
    """Compute features of a window."""

    r1_mean_mapq = 0
    r2_mean_mapq = 0
    r1_total = 0
    r2_total = 0
    yg_mean = 0
    yg_reads = 0
    ya_mean = 0
    ya_reads = 0
    ys_mean = 0
    ys_reads = 0

    for r in window.reads:
        if r.is_read1:
            r1_mean_mapq += r.mapping_quality
            r1_total += 1
        else:
            r2_mean_mapq += r.mapping_quality
            r2_total += 1

        if r.has_tag("YG"):
            yg_mean += r.get_tag("YG")
            yg_reads += 1
        if r.has_tag("YA"):
            ya_mean += r.get_tag("YA")
            ya_reads += 1
        if r.has_tag("YS"):
            ys_mean += r.get_tag("YS")
            ys_reads += 1

    if r1_total:
        r1_mean_mapq /= r1_total
    if r2_total:
        r2_mean_mapq /= r2_total
    if yg_reads:
        yg_mean /= yg_reads
    if ya_reads:
        ya_mean /= ya_reads
    if ys_reads:
        ys_mean /= ys_reads

    return {
        "chr": window.chr,
        "start": window.start,
        "end": window.end,
        "total_reads": r1_total + r2_total,
        "r1_total": r1_total,
        "r2_total": r2_total,
        "r1_mean_mapq": r1_mean_mapq,
        "r2_mean_mapq": r2_mean_mapq,
        "yg_mean": yg_mean,
        "ya_mean": ya_mean,
        "ys_mean": ys_mean,
    }


def window_filter(window):
    """Only yield windows that pass filter."""

    if window["total_reads"] >= 3 and window["ya_mean"] > window["yg_mean"]:
        return True
    else:
        return False


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    # open bam file
    bam = pysam.AlignmentFile(snakemake.input["bam"], "rb")

    # make windows and get features
    df = []

    for contig in bam.references:

        reflen = bam.get_reference_length(contig)

        # make generator for background windows (10kb)
        bg_windows = windows(bam, contig, window_size=10000, step_size=25)

        # get first window
        try:
            bg = next(bg_windows)
            bg_center = bg.start + (bg.end - bg.start) / 2
        except StopIteration:
            continue

        # iterate over windows
        for w in windows(bam, contig, window_size=750, step_size=250):
            f = window_features(w)
            if not window_filter(f):
                continue

            w_center = w.start + (w.end - w.start) / 2
            while True:
                if (
                    (w_center == bg_center)
                    or (w_center < bg_center and bg.start == 0)
                    or (w_center > bg_center and bg.end == reflen)
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

            assert (
                "chr_bg" in f
            ), "failed to find background window for window at {}:{}-{}".format(
                w.chr, w.start, w.end
            )

            df.append(f)

    df = pd.DataFrame(df)

    # add cell_id and donor_id columns
    df["cell_id"] = snakemake.wildcards.sample
    df["donor_id"] = snakemake.wildcards.donor

    # save
    df.to_parquet(snakemake.output[0])

    sys.stderr.close()
