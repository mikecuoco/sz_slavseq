#!/usr/bin/env python
__author__ = "Rohini Gadde"

import pandas as pd
import sys
import pysam

def main():
    df = pd.read_csv(sep="\t", names=["chr", "start", "end"])

    df_label = pd.DataFrame(columns=["chr", "start", "end"])
    samfile = pysam.AlignmentFile(snakemake.input.bam, "rb")

    # for each non-ref insertion...
    for i in df.index:
        chrom = df["chr"][i]
        start = df["start"][i]
        end = df["end"][i]

        # search bulk alignment within +/- 200bp frame of insertion
        bulk = pysam.fetch(chrom, start - 200, end + 200)
        read2 = False

        if len(bulk) > 0:
            for read in bulk:
                read2 = read.is_read2
        
        # keep insertions corresponding to reads in bulk alignment; else, discard
        if read2:
            df_label.loc[len(df_label.index)] = [chrom, start, end]

if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")
    main()
    sys.stderr.close()
