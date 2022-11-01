#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import pysam
import pandas as pd
import sys
from src.genome.Genome import Genome
from src.features.TabixSam import TabixSam
from src.features.ttaaaa import ENSearch
from src.features.WindowFeatures import WindowFeatures
from src.features.occupied_windows import occupied_windows_in_genome


def flank_features(df):
    for i in range(1, 8):
        flank_size = 2**i
        field_name = "flank_" + str(flank_size) + "_max_reads"
        df[field_name] = (
            df["all_reads.count"]
            .rolling(window=2 * flank_size + 1, center=True, min_periods=1)
            .max()
            .fillna(0)
        )


def main(
    filename: str,
    genome: Genome,
    window_size: int,
    window_step: int,
    min_mapq,
    min_ya,
    max_yg,
    min_secondary_mapq,
    library_3_or_5: int,
    ensearch: ENSearch,
):
    # load the alignment
    tbx = TabixSam(pysam.TabixFile(filename))

    # create the windows
    windows = occupied_windows_in_genome(genome, window_size, window_step, filename)

    # first pass through the genome to get features for each window
    df = pd.DataFrame()
    for w in windows:
        wf = wf = WindowFeatures(
            tbx,
            w.chrom,
            w.start,
            w.end,
            min_mapq,
            min_ya,
            max_yg,
            min_secondary_mapq,
            library_3_or_5,
            ensearch,
        )
        f = wf.features()
        f["chrom"] = str(w.chrom)
        f["start"] = w.start
        f["end"] = w.end

        df = df.append(f, ignore_index=True)

    # set index
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df.set_index(["chrom", "start", "end"], inplace=True)

    # second pass through the genome to get collect flanking features
    for i in range(1, 8):
        flank_size = 2**i
        field_name = "flank_" + str(flank_size) + "_max_reads"
        df[field_name] = (
            df["all_reads.count"]
            .rolling(window=2 * flank_size + 1, center=True, min_periods=1)
            .max()
            .fillna(0)
        )

    return df


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    genome = Genome(snakemake.input.chromsizes)
    df = main(
        snakemake.input.bgz,
        genome,
        snakemake.params.window_size,
        snakemake.params.window_step,
        snakemake.params.min_mapq,
        snakemake.params.min_ya,
        snakemake.params.max_yg,
        snakemake.params.min_secondary_mapq,
        snakemake.params.library_3_or_5,
        ENSearch(genome, snakemake.input.fa, 50, 15),
    )

    # save
    df.to_pickle(snakemake.output[0])