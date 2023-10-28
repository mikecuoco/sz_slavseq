#!/usr/bin/env python
# Created on: May 12, 2023 at 2:22:04 PM
__author__ = "Michael Cuoco"

import pysam, sys

sys.stderr = open(snakemake.log[0], "w")

pysam.set_verbosity(0)

# initialize variables
r1, r2 = type("", (), {})(), type("", (), {})()
r1.query_name, r2.query_name = "", None
assert r1.query_name != r2.query_name

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
                while (r1.query_name != r2.query_name) or (
                    r1.query_name != r_l1.query_name
                ):
                    r_genome = next(genome_bam)

                    # get each read's primary alignment
                    if r_genome.is_secondary or r_genome.is_supplementary:
                        continue

                    # get the read, save as r1 or r2
                    if r_genome.is_read1:
                        r1 = r_genome
                    else:
                        r2 = r_genome

                # skip if read1 is unmapped
                if r1.is_unmapped:
                    continue

                assert (
                    r1.query_name == r2.query_name == r_l1.query_name
                ), "Read IDs do not match"

                # add tags to read1
                r1.set_tag("L1", r_l1.get_tag("AS"))  # L1 alignment score
                r1.set_tag("MS", r2.get_tag("AS"))  # mate alignment score
                r1.set_tag("LS", r_l1.reference_start)  # L1 start
                r1.set_tag("LE", r_l1.reference_end)  # L1 end
                r1.set_tag("LA", r_l1.query_sequence.count("A"))  # L1 A count
                r1.set_tag(
                    "ML",
                    r2.infer_read_length()
                    if r2.infer_read_length()
                    else len(r2.query_sequence),
                )  # mate read length

                # write reads to output bam
                out_bam.write(r1)
                out_bam.write(r2)

sys.stderr.close()
