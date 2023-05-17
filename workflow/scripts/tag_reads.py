#!/usr/bin/env python
# Created on: May 12, 2023 at 2:22:04 PM
__author__ = "Michael Cuoco"

import pysam, sys

sys.stderr = open(snakemake.log[0], "w")

pysam.set_verbosity(0)

# initialize variables
r1, r2 = None, None

# iterate over read IDs in line1_bam, match with genome_bam reads
# save primary r1 from genome, primary r2 from genome, and primary r2 from line-1 to tab-delimited file
with pysam.AlignmentFile(snakemake.input["genome_bam"], "rb") as genome_bam:
    with pysam.AlignmentFile(snakemake.input["line1_bam"], "rb") as line1_bam:
        with pysam.AlignmentFile(
            snakemake.output[0], "wb", template=genome_bam
        ) as out_bam:
            # iterate over alignments in line1_bam
            for r_l1 in line1_bam:
                # skip secondary and supplementary alignments
                if r_l1.is_secondary or r_l1.is_supplementary:
                    continue

                # find matching reads in genome_bam
                while True:
                    r_genome = next(genome_bam)

                    # skip secondary and supplementary alignments
                    if r_genome.is_secondary or r_genome.is_supplementary:
                        continue

                    # skip if genome read ID does not match line1 read ID
                    if r_genome.query_name != r_l1.query_name:
                        continue

                    # get the read, save as r1 or r2
                    if r_genome.is_read1:
                        r1 = r_genome
                    else:
                        r2 = r_genome

                    # skip if read1 or read2 is None
                    if r1 is None or r2 is None:
                        continue

                    # break loop when both read1 and read2 are found
                    elif r1.query_name == r2.query_name:
                        break

                # skip if read1 is a duplicate
                if r1.is_duplicate or (r1.is_unmapped and r2.is_unmapped):
                    continue

                # add tags to read1
                r1.set_tag("ML", r_l1.get_tag("AS"))
                r1.set_tag("MG", r2.get_tag("AS"))
                r1.set_tag("MS", r_l1.reference_start)
                r1.set_tag("MA", r_l1.query_sequence.count("A"))

                # write read1 to output bam
                out_bam.write(r1)

sys.stderr.close()
