#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import functools, sys, os
import polars as pl
import pandas as pd

# get access to the src directory
sys.path.append((os.path.abspath("workflow")))
from src.genome.Genome import Genome
from src.genome import interval_generator as ig
from src.features.occupied_windows import occupied_windows_in_genome
from src.genome.Interval import Interval

os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)


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
        "in_rmsk",
        snakemake.params.window_size,
        snakemake.params.window_step,
    )
    return pl.DataFrame(df.reset_index())


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
    return pl.DataFrame(df.reset_index())


def label(x):
    if x["in_rmsk"]:
        return "RL1"
    elif x["in_NRdb"]:
        return "KNRGL"
    else:
        return "OTHER"


def get_germline_l1(
    tbx_fn: str,
    genome: Genome,
    window_size: int,
    window_step: int,
):

    # create the windows
    windows = occupied_windows_in_genome(genome, window_size, window_step, tbx_fn)

    w_list = []
    for w in windows:
        if w is None:
            continue
        w_list.append(w.as_tuple())
    df = pl.DataFrame(w_list, schema=["chrom", "start", "end"])

    # merge ref and nonref l1
    df = (
        df.join(read_non_ref_db(), on=["chrom", "start", "end"], how="left")
        .join(read_reference_l1(), on=["chrom", "start", "end"], how="left")
        .fill_null(False)
    )

    # only keep windows with L1 in either the reference or nonref database
    df = df.filter(pl.col("in_NRdb") | pl.col("in_rmsk"))

    # if ref and nonref l1 are present at the same site, set nonref to false
    df = df.with_column(
        pl.when(pl.col("in_rmsk") & pl.col("in_NRdb"))
        .then(False)
        .otherwise(pl.col("in_NRdb"))
        .alias("in_NRdb")
    )

    return df


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    # get germline L1 from bulk data
    germline_df = get_germline_l1(
        snakemake.input.bgz[0],
        Genome(snakemake.input.chromsizes),
        snakemake.params.window_size,
        snakemake.params.window_step,
    )

    # save
    germline_df.with_columns(
        [
            pl.struct(pl.col(["in_NRdb", "in_rmsk"])).apply(label).alias("label"),
            pl.lit(snakemake.wildcards.db).alias("db"),
        ]
    ).drop(["in_NRdb", "in_rmsk"]).write_csv(snakemake.output.bulk)

    # get features from single cell data
    features_df = pl.concat([pl.read_parquet(f) for f in snakemake.input.features])

    # merge
    df = features_df.join(
        germline_df, on=["chrom", "start", "end"], how="left"
    ).fill_null(False)

    # collapse to single column
    df = df.with_column(
        pl.struct(pl.col(["in_NRdb", "in_rmsk"])).apply(label).alias("label")
    ).drop(["in_NRdb", "in_rmsk"])

    if snakemake.params.rmL1:
        df = df.filter(pl.col("label") != "RL1")

    # add db and build columns
    df = df.with_column(pl.lit(snakemake.wildcards.db).alias("db"))

    # save
    df.write_parquet(snakemake.output.mda)
    sys.stderr.close()
