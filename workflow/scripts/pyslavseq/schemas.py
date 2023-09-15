import numpy as np
from collections import namedtuple

# for converting read to a hashable format
Read = namedtuple(
    "Read",
    [
        "reference_name",
        "reference_start",
        "reference_end",
        "read_length",
        "is_read1",
        "is_read2",
        "is_forward",
        "is_reverse",
        "mapping_quality",
        "num_supp_alignments",
        "cigar_string",
        "alignment_score",
        "L1_alignment_score",
        "L1_reference_start",
        "L1_reference_end",
        "L1_Acount",
        "mate_alignment_score",
        "mate_read_length",
        "mate_mapping_quality",
        "mate_cigar_string",
        "mate_reference_name",
        "mate_reference_start",
        "is_proper_pair",
        "isref_read",
    ],
)

# define schema for feature tables
TAGS = [
    "alignment_score",
    "alignment_score_normed",
    "L1_alignment_score",
    "L1_alignment_score_normed",
    "L1_reference_start",
    "L1_reference_end",
    "L1_Acount",
    "mate_alignment_score",
    "mate_alignment_score_normed",
    "mate_read_length",
    "num_supp_alignments",
]

# define schema for peaks
PEAKS_SCHEMA = {
    "Chromosome": "",
    "Start": np.int32(),
    "End": np.int32(),
    "n_ref_reads": np.int32(),
    "n_reads": np.int32(),
    "rpm": np.float32(),
    "max_mapq": np.int32(),
    "n_unique_starts": np.int32(),
}
