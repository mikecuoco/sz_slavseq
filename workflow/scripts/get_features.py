#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import pysam
import pandas as pd
import sys
from statistics import mean
from src.genome.Genome import Genome
from src.genome import interval_generator as ig

def mean_frag_len(reads):
    l = [r.template_length for r in reads if r.is_paired and r.is_read1]
    return 0 if len(l) == 0 else abs(mean(l))


def mean_r2_l1_score(reads):
    l = [r.get_tag("YA") for r in reads if r.is_read2 and r.has_tag("YA")]
    return 0 if len(l) == 0 else mean(l)


def mean_r2_ref_score(reads):
    l = [r.get_tag("YG") for r in reads if r.is_read2 and r.has_tag("YG")]
    return 0 if len(l) == 0 else mean(l)


def n_r2_uniq_starts(reads):
    return len(set([r.reference_start for r in reads if r.is_read2]))


def n_r1_uniq_starts(reads):
    return len(set([r.reference_start for r in reads if r.is_read1]))


def cell_features(
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
        # get reads in window
        reads = [r for r in bam.fetch(w.chrom, w.start, w.end)]
        if len(reads) < min_reads:
            continue
        f = {}
        f["chrom"] = str(w.chrom)
        f["start"] = w.start
        f["end"] = w.end
        f["count"] = len(reads)

        f["mean_frag_len"] = mean_frag_len(reads)
        f["n_r1_uniq_starts"] = n_r1_uniq_starts(reads)
        f["n_r2_uniq_starts"] = n_r2_uniq_starts(reads)
        f["mean_r2_l1_score"] = mean_r2_l1_score(reads)
        f["mean_r2_ref_score"] = mean_r2_ref_score(reads)

        w_list.append(f)

    df = pd.DataFrame(w_list)
    # set index
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df.set_index(["chrom", "start", "end"], inplace=True)

    return df


def cell_flank_features(df):

    df.sort_index(inplace=True)

    for i in range(1, 8):
        flank_size = 2**i
        field_name = "flank_" + str(flank_size) + "_max_reads"
        df[field_name] = (
            df["count"]
            .rolling(window=2 * flank_size + 1, center=True, min_periods=1)
            .max()
            .fillna(0)
        )

    return df


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    # first pass through the genome to get features for each window
    df = cell_features(
        snakemake.input.bam,
        Genome(snakemake.input.chromsizes),
        snakemake.params.window_size,
        snakemake.params.window_step,
        snakemake.params.min_reads,
    )

    # second pass through the genome to get collect flanking features
    df = cell_flank_features(df)

    # write the features to a file
    df["cell_id"] = snakemake.wildcards.sample
    df["donor_id"] = snakemake.wildcards.donor
    df.set_index(["cell_id", "donor_id"], append=True, inplace=True)
    df.to_pickle(snakemake.output[0])

    sys.stderr.close()
