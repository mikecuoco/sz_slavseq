#!/usr/bin/env python
# Created on: Nov 21, 2023 at 12:13:28 PM
__author__ = "Michael Cuoco"
# DEPECATED!!

import logging

logger = logging.getLogger(__name__)  # configure logging

from time import time
from collections import namedtuple
from pysam import AlignedSegment
import numpy as np


# use named tuple to convert AlignedSegment
# gives a speedup over using pysam.AlignedSegment (about 1.5x)
Read = namedtuple(
    "Read",
    [
        "query_name",
        "reference_name",
        "reference_start",
        "reference_end",
        "five_end",
        "three_end",
        "breakpoint_end",
        "read_length",
        "is_read1",
        "is_read2",
        "is_forward",
        "is_reverse",
        "is_secondary",
        "is_supplementary",
        "is_mapped",
        "mapping_quality",
        "num_supp_alignments",
        "cigarstring",
        "three_end_clippedA",
        "three_end_clipped_length",
        "alignment_score",
        "L1_alignment_score",
        "L1_mapping_quality",
        "L1_reference_start",
        "L1_reference_end",
        "L1_Acount",
        "mate_alignment_score",
        "mate_read_length",
        "mate_mapping_quality",
        "mate_cigarstring",
        "mate_reference_name",
        "mate_reference_start",
        "is_proper_pair",
        "is_ref_read",
        "is_duplicate",
    ],
)


def read_to_namedtuple(read: AlignedSegment) -> Read:
    "convert pysam.AlignedSegment to hashable namedtuple Read"
    return Read(
        read.query_name,
        read.reference_name,
        read.reference_start,
        read.reference_end,
        read.reference_end if read.is_reverse else read.reference_start,
        read.reference_start if read.is_reverse else read.reference_end,
        read.reference_end if read.is_read2 ^ read.is_reverse else read.reference_start,
        read.infer_read_length(),
        read.is_read1,
        read.is_read2,
        not read.is_reverse,
        read.is_reverse,
        read.is_secondary,
        read.is_supplementary,
        read.is_mapped,
        read.mapping_quality,
        len(read.get_tag("SA").split(";")[:-1]) if read.has_tag("SA") else 0,  # type: ignore
        read.cigarstring,
        read.get_tag("CA") if read.has_tag("CA") else None,
        read.get_tag("CL") if read.has_tag("CL") else None,
        read.get_tag("AS") if read.has_tag("AS") else None,
        read.get_tag("L1") if read.has_tag("L1") else None,
        read.get_tag("LQ") if read.has_tag("LQ") else None,
        read.get_tag("LS") if read.has_tag("LS") else None,
        read.get_tag("LE") if read.has_tag("LE") else None,
        read.get_tag("LA") if read.has_tag("LA") else None,
        read.get_tag("MS") if read.has_tag("MS") else None,
        read.get_tag("ML") if read.has_tag("ML") else None,
        read.get_tag("MQ") if read.has_tag("MQ") else None,
        read.get_tag("MC") if read.has_tag("MC") else None,
        read.next_reference_name,
        read.next_reference_start,
        read.is_proper_pair,
        read.get_tag("RR"),
        read.is_duplicate,
    )
