#!/usr/bin/env python
# adapted from Ricardo's jupyter notebook from rf pipeline: pipeline/eul1db/to_windows/_h/main.ipynb
__author__ = "Michael Cuoco, Rohini Gadde"

import pandas as pd
from .Genome import Genome
from .Interval import Interval
from . import interval_generator as ig
import sys, gc, traceback


def read_rmsk(rmsk_outfile):
    """
    Read the repeatmasker output table and return locations of L1HS and L1PA2-6
    """
    # read the rmsk file
    df0 = pd.read_csv(
        rmsk_outfile,
        skiprows=3,
        delim_whitespace=True,
        names=["chrom", "start", "end", "strand", "repeat"],
        usecols=[4, 5, 6, 8, 9],
    )

    # filter for rep_names
    rep_names = ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]
    # logging.info(f"Filtering for rep_names: {rep_names}")
    df0 = df0[df0["repeat"].isin(rep_names)]

    # save to new dataframe
    df1 = pd.DataFrame()
    df1["chrom"] = df0["chrom"]
    # set start positions depending on strand
    df1["start"] = df0.apply(
        lambda x: x["end"] if x["strand"] != "+" else x["start"], axis=1
    )
    df1["end"] = df1["start"]
    df1["start"] -= 1  # make zero-based

    return df1


def make_l1_windows(df, chromsizes, field):
    l1_pos = set()

    for (_, chrom, start, end) in df[["chrom", "start", "end"]].itertuples():
        l1_pos.update([Interval(chrom, start, end)])

    genome = Genome(chromsizes)
    xx = list(ig.windows_overlapping_intervals(genome, l1_pos, 750, 250))

    l1_df = pd.DataFrame.from_records(
        (x.as_tuple() for x in xx), columns=["chrom", "start", "end"]
    ).set_index(["chrom", "start", "end"])
    l1_df[field] = True

    return l1_df


def make_genome_windows(chromsizes):
    genome = Genome(chromsizes)
    wg_iter = ig.windows_in_genome(genome, 750, 250)
    df = pd.DataFrame.from_records(
        (xx.as_tuple() for xx in wg_iter), columns=["chrom", "start", "end"]
    ).set_index(["chrom", "start", "end"])
    return df
