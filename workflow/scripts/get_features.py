#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import pysam
import polars as pl
import sys, os

# get access to the src directory
sys.path.append((os.path.abspath("workflow")))
from src.genome.Genome import Genome
from src.features.TabixSam import TabixSam
from src.features.WindowFeatures import WindowFeatures
from src.features.occupied_windows import occupied_windows_in_genome

os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)


def flank_features(df):

    # sort
    df = df.sort(["chrom", "start", "end"]).to_pandas()

    for i in range(1, 8):
        flank_size = 2**i
        field_name = "flank_" + str(flank_size) + "_max_reads"
        df[field_name] = (
            df["all_reads.count"]
            .rolling(window=2 * flank_size + 1, center=True, min_periods=1)
            .max()
            .fillna(0)
        )

    return pl.DataFrame(df)


def features(
    filename: str,
    genome: Genome,
    window_size: int,
    window_step: int,
    min_mapq,
    min_ya,
    max_yg,
    min_secondary_mapq,
    library_3_or_5: int,
    ensearch,
):
    # load the alignment
    tbx = TabixSam(pysam.TabixFile(filename))

    # create the windows
    windows = occupied_windows_in_genome(genome, window_size, window_step, filename)

    # first pass through the genome to get features for each window
    w_list = []
    for w in windows:
        wf = WindowFeatures(
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
        w_list.append(f)

    df = pl.DataFrame(w_list)

    return df


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    df = features(
        snakemake.input.bgz,
        Genome(snakemake.input.chromsizes),
        snakemake.params.window_size,
        snakemake.params.window_step,
        snakemake.params.min_mapq,
        snakemake.params.min_ya,
        snakemake.params.max_yg,
        snakemake.params.min_secondary_mapq,
        snakemake.params.library_3_or_5,
        None,
    )

    df = flank_features(df)

    # add cell_id and donor_id columns
    df = df.with_column(pl.lit(snakemake.wildcards.sample).alias("cell_id"))
    df = df.with_column(pl.lit(snakemake.wildcards.donor).alias("donor_id"))

    # save
    df.write_parquet(snakemake.output[0])

    sys.stderr.close()
