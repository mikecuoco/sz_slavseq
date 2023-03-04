#!/usr/bin/env python3
"""This program takes headerless paired-end SAM input sorted by read name. For a read pair, the main alignment is defined as the primary alignment of R1 (i.e., the alignment not having flags 2048 or 256). The coordinates of that alignment are used to index all other alignments belonging to the pair (i.e. R2 and secondary alignments) using a tabix index.

Example:
id1/r1: maps to chr1:1000-1080 (primary r1 alignment)
id1/r1: maps to chr14:2020-2040 (secondary r1 alignment)
id1/r2: maps to chrX:10-110 (primary r2 alignment)

This way, a search for chr1:1000-2000 will return all 3 alignments. A
search for chr14:2020-2040 or chrX:10-110 will return nothing.
"""

import sys, pysam

def output_alignments(alignments, primary_r1, qname):
    if len(alignments) > 0 and primary_r1:
        for al in alignments:
            joined = "\t".join(
                [
                    primary_r1.reference_name,
                    str(primary_r1.reference_start),
                    str(primary_r1.reference_end),
                    al.to_string(),
                ]
            )
            print(joined, file=sys.stdout)


if __name__ == "__main__":
    if len(sys.argv) > 2:
        print("Usage: sam_to_tabix.py input.bam > output.tabix", file=sys.stderr)
        sys.exit(1)


    # iterate over all alignments
    with pysam.AlignmentFile(str(sys.argv[1]), "rb") as bam:

        # initialize 
        alignments = []
        b = bam.fetch(until_eof=True)
        primary_r1 = None
        for al in b:
            if al.is_unmapped:
                continue
            last_qname = al.query_name
            alignments.append(al)
            if al.is_read1 and (not al.is_secondary) and (not al.is_supplementary):
                primary_r1 = al
                break

        # iterate over all alignments
        kept, skipped = 0,0
        for al in b:
            if al.is_unmapped:
                continue

            # grab the read name
            current_qname = al.query_name

            # check if read has new name
            if last_qname != current_qname:
                output_alignments(alignments, primary_r1, last_qname)
                if primary_r1:
                    kept += len(alignments) 
                else:
                    skipped += len(alignments)
                alignments = []
                primary_r1 = None
                last_qname = current_qname

            # add the alignment to the list
            alignments.append(al)

            # check if this is the primary r1 alignment, if so, store it
            if al.is_read1 and (not al.is_secondary) and (not al.is_supplementary):
                assert primary_r1 is None, "Primary R1 error"
                primary_r1 = al

output_alignments(alignments, primary_r1, last_qname)
if primary_r1:
    kept += len(alignments) 
else:
    skipped += len(alignments)
print(kept, "alignments processed", file=sys.stderr)
print(skipped, "alignments without primary read1 were skipped", file=sys.stderr)
