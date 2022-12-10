#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import functools, sys
import pandas as pd
import pysam
from src.genome.Genome import Genome
from src.genome import interval_generator as ig
from src.genome.windows import read_rmsk, make_l1_windows


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


def get_germline_l1(
    bam_fn: str,
    genome: Genome,
    window_size: int,
    window_step: int,
    min_reads: int,
):
    # create the windows
    windows = ig.windows_in_genome(genome, window_size, window_step)
    bam = pysam.AlignmentFile(bam_fn, "rb")

    # iterate over the windows to compute the features
    w_list = []
    for w in windows:
        reads = len([r for r in bam.fetch(w.chrom, w.start, w.end)])
        if reads >= min_reads:
            w_list.append(w.as_tuple())
    df = pd.DataFrame(w_list, columns=["chrom", "start", "end"])

    # set index
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df.set_index(["chrom", "start", "end"], inplace=True)
    df = (
        df.merge(read_non_ref_db(), left_index=True, right_index=True, how="left")
        .merge(read_reference_l1(), left_index=True, right_index=True, how="left")
        .fillna({"in_NRdb": False, "reference_l1hs_l1pa2_6": False})
    )
    df = df.loc[df["in_NRdb"] | df["reference_l1hs_l1pa2_6"]]
    df.loc[
        df["in_NRdb"] & df["reference_l1hs_l1pa2_6"], ["in_NRdb"]
    ] = False  # if ref and nonref l1 are present at the same site, set nonref to false
    return df


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    # get germline L1 from bulk data
    germline_df = get_germline_l1(
        snakemake.input.bam[0],
        Genome(snakemake.input.chromsizes),
        snakemake.params.window_size,
        snakemake.params.window_step,
        snakemake.params.min_reads,
    )

    # get features from single cell data
    features_df = pd.concat([pd.read_pickle(f) for f in snakemake.input.features])

    # merge and return
    df = features_df.merge(
        germline_df, left_index=True, right_index=True, how="left"
    ).fillna({"in_NRdb": False, "reference_l1hs_l1pa2_6": False})
    df.to_pickle(snakemake.output[0])
    sys.stderr.close()
