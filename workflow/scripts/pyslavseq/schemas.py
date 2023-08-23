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

FEATURES_SCHEMA = {
    "Chromosome": "",
    "Start": np.int32(),
    "End": np.int32(),
    "n_fwd": np.int32(),
    "n_rev": np.int32(),
    "n_proper_pairs": np.int32(),
    "n_ref_reads": np.int32(),
    "3end_gini": np.float32(),
    "5end_gini": np.float32(),
    "max_mapq": np.int8(),
    "n_reads": np.int32(),
    "rpm": np.float64(),
    "orientation_bias": np.float32(),
    "frac_proper_pairs": np.float32(),
}

for tag in TAGS:
    for q in [0, 0.25, 0.5, 0.75, 1]:
        FEATURES_SCHEMA[f"{tag}_q{q}"] = np.float32()
    FEATURES_SCHEMA[f"{tag}_mean"] = np.float32()


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
