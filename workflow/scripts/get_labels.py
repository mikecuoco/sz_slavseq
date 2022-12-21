#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import functools, sys
import pandas as pd
import pysam
from src.genome.Genome import Genome
from src.features.TabixSam import TabixSam
from src.genome import interval_generator as ig
from src.features.occupied_windows import occupied_windows_in_genome
from src.genome.Interval import Interval

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
    rep_names = [
        "L1HS_3end",
        "L1PA2_3end",
        "L1PA3_3end",
        "L1PA4_3end",
        "L1PA5_3end",
        "L1PA6_3end",
    ]
    # logging.info(f"Filtering for rep_names: {rep_names}")
    df0 = df0[df0["repeat"].isin(rep_names)]

    # save to new dataframe
    df1 = pd.DataFrame()
    df1["chrom"] = df0["chrom"].astype(str)
    # set start positions depending on strand
    df1["start"] = df0.apply(
        lambda x: x["end"] if x["strand"] != "+" else x["start"], axis=1
    )
    df1["end"] = df1["start"]
    df1["start"] -= 1  # make zero-based

    return df1


def make_l1_windows(df, chromsizes, field, window_size, window_step):
    l1_pos = set()

    for (_, chrom, start, end) in df[["chrom", "start", "end"]].itertuples():
        l1_pos.update([Interval(chrom, start, end)])

    genome = Genome(chromsizes)
    xx = list(
        ig.windows_overlapping_intervals(genome, l1_pos, window_size, window_step)
    )

    l1_df = pd.DataFrame.from_records(
        (x.as_tuple() for x in xx), columns=["chrom", "start", "end"]
    ).set_index(["chrom", "start", "end"])
    l1_df[field] = True

    return l1_df


@functools.lru_cache()
def read_reference_l1():
    df = read_rmsk(snakemake.input.ref_l1)
    df = make_l1_windows(
        df,
        snakemake.input.chromsizes,
        "reference_l1hs_l1pa2_6",
        snakemake.params.window_size,
        snakemake.params.window_step,
    )
    return df


@functools.lru_cache()
def read_non_ref_db():
    df = pd.read_csv(
        snakemake.input.non_ref_l1,
        sep="\t",
        header=None,
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": int, "end": int},
    )
    df = make_l1_windows(
        df,
        snakemake.input.chromsizes,
        "in_NRdb",
        snakemake.params.window_size,
        snakemake.params.window_step,
    )
    return df


def label(df):
    for x in df[["in_NRdb", "reference_l1hs_l1pa2_6"]].itertuples():
        if x.reference_l1hs_l1pa2_6:
            yield "RL1"
        elif x.in_NRdb:
            yield "KNRGL"
        else:
            yield "OTHER"


def get_germline_l1(
    tbx_fn: str,
    genome: Genome,
    window_size: int,
    window_step: int,
    min_reads: int,
    db: str,
):

    tbx = TabixSam(pysam.TabixFile(tbx_fn))

    # create the windows
    windows = occupied_windows_in_genome(genome, window_size, window_step, tbx_fn)

    w_list = []
    for w in windows:
        reads = len([r for r in tbx._fetch(w.chrom, w.start, w.end)])
        if reads >= min_reads:
            w_list.append(w.as_tuple())
    df = pd.DataFrame(w_list, columns=["chrom", "start", "end"])

    # set index
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df.set_index(["chrom", "start", "end"], inplace=True)

    # merge ref and nonref l1
    df = (
        df.merge(read_non_ref_db(), left_index=True, right_index=True, how="left")
        .merge(read_reference_l1(), left_index=True, right_index=True, how="left")
        .fillna({"in_NRdb": False, "reference_l1hs_l1pa2_6": False})
    )

    df = df.loc[df["in_NRdb"] | df["reference_l1hs_l1pa2_6"]]
    df.loc[
        df["in_NRdb"] & df["reference_l1hs_l1pa2_6"], ["in_NRdb"]
    ] = False  # if ref and nonref l1 are present at the same site, set nonref to false

    # add database source
    df["db"] = df["in_NRdb"].apply(lambda x: db if x == True else "rmsk")

    return df


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    # get germline L1 from bulk data
    germline_df = get_germline_l1(
        snakemake.input.bgz[0],
        Genome(snakemake.input.chromsizes),
        snakemake.params.window_size,
        snakemake.params.window_step,
        snakemake.params.min_reads,
        snakemake.wildcards.db,
    )

    # get features from single cell data
    features_df = pd.concat([pd.read_pickle(f) for f in snakemake.input.features])

    # merge and return
    df = features_df.merge(
        germline_df, left_index=True, right_index=True, how="left"
    ).fillna({"in_NRdb": False, "reference_l1hs_l1pa2_6": False, "db": False})

    # collapse to single column
    df["label"] = pd.Series(label(df), index=df.index)
    df.drop(["in_NRdb", "reference_l1hs_l1pa2_6"], axis=1, inplace=True)

    df["build"] = snakemake.wildcards.ref

    # error check
    assert len(set(df["label"])) == 3, "Not all labels are present"

    # save
    df.to_pickle(snakemake.output[0])
    sys.stderr.close()
