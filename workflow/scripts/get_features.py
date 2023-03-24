#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import pysam
import pandas as pd
from collections import deque, namedtuple, defaultdict
import sys


def reads(bam, contig):
    """Yield reads from a bam file for a given contig."""
    for r in bam.fetch(contig, multiple_iterators=True):
        if not r.is_duplicate and not r.is_unmapped and r.mapping_quality > 0:
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

    f = defaultdict(int)

    for r in window.reads:
        if r.is_read1:
            f["r1_mean_mapq"] += r.mapping_quality
            if r.is_reverse:
                f["r1_rev"] += 1
            else:
                f["r1_fwd"] += 1
        else:
            f["r2_mean_mapq"] += r.mapping_quality
            if r.is_reverse:
                f["r2_rev"] += 1
            else:
                f["r2_fwd"] += 1

        if r.has_tag("YG"):
            f["yg_mean"] += r.get_tag("YG")
            f["yg_reads"] += 1
        if r.has_tag("YA"):
            f["ya_mean"] += r.get_tag("YA")
            f["ya_reads"] += 1
        if r.has_tag("YS"):
            f["ys_mean"] += r.get_tag("YS")
            f["ys_reads"] += 1

    f["r1_total"] = f["r1_rev"] + f["r1_fwd"]
    if f["r1_total"] > 0:
        f["r1_mean_mapq"] /= f["r1_total"]
    f["r2_total"] = f["r2_rev"] + f["r2_fwd"]
    if f["r2_total"] > 0:
        f["r2_mean_mapq"] /= f["r2_total"]
    if f["yg_reads"]:
        f["yg_mean"] /= f["yg_reads"]
    if f["ya_reads"]:
        f["ya_mean"] /= f["ya_reads"]
    if f["ys_reads"]:
        f["ys_mean"] /= f["ys_reads"]

    f["chr"] = window.chr
    f["start"] = window.start
    f["end"] = window.end
    f["total_reads"] = f["r1_total"] + f["r2_total"]

    return f


def window_filter(window):
    """Only yield windows that pass filter."""

    # if window["total_reads"] >= 3 and window["ya_mean"] > window["yg_mean"]:
    if window["total_reads"] >= 3:
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
