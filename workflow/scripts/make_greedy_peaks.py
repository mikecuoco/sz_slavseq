#!/usr/bin/env python
# Created on: Jul 10, 2024 at 10:36:44 PM
__author__ = "Michael Cuoco"


import logging

logger = logging.getLogger(__name__)

import pyBigWig
import pyranges as pr
from collections import defaultdict
import pandas as pd


def greedy(windows, blacklist_size: int = 1000):
    """
    Greedily choose windows with maximum number of reads
    :param windows: generator of windows from a bigwig file of a chromosome
    :param blacklist_size: size of region to exclude around chosen window
    """

    # exhaust windows generator and sort by number of deduplicated reads
    windows = sorted(windows, key=lambda w: w["n_reads"], reverse=True)
    blacklist = []  # initialize list of regions to exclude
    for w in windows:
        # if the window is not in the blacklist, yield it
        if not any(((w["Start"] <= x[1]) and (w["End"] >= x[0])) for x in blacklist):
            yield w
            blacklist.append((w["Start"] - blacklist_size, w["End"] + blacklist_size))

    logger.info(
        f"Found {len(blacklist)} windows using {blacklist_size} bp flanking exclusion zones."
    )


def make_greedy_peaks(bw: str, blacklist_size: int, extend_bp: int) -> pd.DataFrame:
    """
    Make greedy peaks from a bigwig file
    :param bw: path to bigwig file
    :param blacklist_size: size of region to exclude around each max position in greedy algorithm
    :param extend_bp: number of base pairs to extend the peak in each direction
    """

    bw = pyBigWig.open(bw)

    # load entire bigwig into mem
    peaks = defaultdict(list)
    tup_to_dict = lambda x: dict(zip(["Start", "End", "n_reads"], x))
    for c in bw.chroms():
        for w in greedy(
            map(tup_to_dict, bw.intervals(c)), blacklist_size=blacklist_size
        ):
            peaks[c].append(w)
        peaks[c] = sorted(peaks[c], key=lambda x: x["Start"])

    # concantenate results across chromosomes, convert to bigwig
    df = pd.concat({k: pd.DataFrame(v) for k, v in peaks.items()}, axis=0)
    df = df.drop(columns="n_reads")
    df = (
        df.reset_index(level=0)
        .rename(columns={"level_0": "Chromosome"})
        .reset_index(drop=True)
    )

    return pr.PyRanges(df).sort().extend(extend_bp).df


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )
    logger = logging.getLogger(__name__)

    df = make_greedy_peaks(snakemake.input.bw, snakemake.params.blacklist_bp, snakemake.params.extend_bp)  # type: ignore
    df.to_csv(snakemake.output.bed, sep="\t", index=False, header=False)
