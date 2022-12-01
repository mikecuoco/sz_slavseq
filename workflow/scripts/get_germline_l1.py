#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import functools, sys
import pandas as pd
import pysam
from src.genome.Genome import Genome
from src.features.occupied_windows import occupied_windows_in_genome
from src.genome.windows import read_rmsk, make_l1_windows
import pdb


@functools.lru_cache()
def read_reference_l1():
    df = read_rmsk(snakemake.input.ref_l1)
    df = make_l1_windows(df, snakemake.input.chromsizes, "reference_l1hs_l1pa2_6")
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
    df = make_l1_windows(df, snakemake.input.chromsizes, "in_NRdb")
    return df


def main(
    filename: str,
    genome: Genome,
    window_size: int,
    window_step: int,
):

    # create the windows
    # TODO use min-reads here
    windows = occupied_windows_in_genome(genome, window_size, window_step, filename)
    df = pd.DataFrame.from_records(
        [w.as_tuple() for w in windows],
        columns=["chrom", "start", "end"],
    )
    df.chrom = df.chrom.astype(str)
    df.set_index(["chrom", "start", "end"], inplace=True)
    df = (
        df.merge(read_non_ref_db(), left_index=True, right_index=True, how="left")
        .merge(read_reference_l1(), left_index=True, right_index=True, how="left")
        .fillna({"in_NRdb": False, "reference_l1hs_l1pa2_6": False})
    )
    df = df.loc[df["in_NRdb"] | df["reference_l1hs_l1pa2_6"]]
    df.loc[df["in_NRdb"] & df["reference_l1hs_l1pa2_6"], ["in_NRdb"]] = False # if ref and nonref l1 are present at the same site, set nonref to false 
    return df


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    genome = Genome(snakemake.input.chromsizes)
    df = main(
        snakemake.input.bgz[0],
        genome,
        snakemake.params.window_size,
        snakemake.params.window_step,
    )
    df.to_pickle(snakemake.output[0])
    sys.stderr.close()
