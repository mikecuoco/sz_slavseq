#!/usr/bin/env python
# Created on: May 12, 2023 at 2:22:04 PM
__author__ = "Michael Cuoco"

import pysam, sys, re, pdb

sys.stderr = open(snakemake.log[0], "w")

pysam.set_verbosity(0)


def is_ref_read(read: pysam.AlignedSegment) -> bool:
    "return True if read is ref, False if non-ref based on mate tag"

    assert not read.is_read2, "Read cannot be read2"
    # get cigar at start of read, accounting for mate and orientation
    if read.is_read1:
        if not read.has_tag("MC"):
            return False
        cigar = re.findall(r"(\d+)([MIDNSHP=X])", str(read.get_tag("MC")))
        end = cigar[-1] if not read.is_reverse else cigar[0]
        clipped = int(end[0]) if end[1] == "H" or end[1] == "S" else 0
        return read.is_proper_pair and (clipped < 30)
    else:
        cigar = re.findall(r"(\d+)([MIDNSHP=X])", read.cigarstring)
        end = cigar[0] if read.is_reverse else cigar[-1]
        clipped = int(end[0]) if end[1] == "H" or end[1] == "S" else 0
        return clipped < 30


def reset_reads():
    """Return empty r1 and r2"""
    r1, r2 = type("", (), {})(), type("", (), {})()
    r1.query_name, r2.query_name = "", None
    return r1, r2


r1, r2 = reset_reads()


def next_primary(bam):
    """Return the next primary alignment in bam"""
    while True:
        r = next(bam)
        if r.is_secondary or r.is_supplementary:
            continue
        else:
            return r


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
                if (
                    r_l1.is_secondary
                    or r_l1.is_supplementary
                    or r_l1.is_unmapped
                    or r_l1.is_read1
                ):
                    continue

                # find matching reads in genome_bam
                while (r1.query_name != r2.query_name) or (
                    r1.query_name != r_l1.query_name
                ):
                    r_genome = next_primary(genome_bam)

                    # get the read, save as r1 or r2
                    if r_genome.is_read1:
                        r1 = r_genome
                    elif r_genome.is_read2:
                        r2 = r_genome
                    else:
                        r1 = r_genome
                        r2 = r_genome

                # skip if read1 is unmapped
                if r1.is_unmapped:
                    continue

                assert (
                    r1.query_name == r2.query_name == r_l1.query_name
                ), "Read IDs do not match"
                r1.set_tag("RR", int(is_ref_read(r1)))
                r1.set_tag("LQ", r_l1.mapping_quality)  # L1 mapping quality
                r1.set_tag("L1", r_l1.get_tag("AS"))  # L1 alignment score
                r1.set_tag("LS", r_l1.reference_start)  # L1 start
                r1.set_tag("LE", r_l1.reference_end)  # L1 end
                r1.set_tag("LA", r_l1.query_sequence.count("A"))  # L1 A count
                r1.set_tag("MS", r2.get_tag("AS"))  # mate alignment score
                r1.set_tag(
                    "ML",
                    r2.infer_read_length()
                    if r2.infer_read_length()
                    else len(r2.query_sequence),
                )  # mate read length

                # write reads to output bam
                out_bam.write(r1)
                if r1.is_read1:
                    out_bam.write(r2)

sys.stderr.close()
