#!/usr/bin/env python
# Created on: Nov 21, 2023 at 12:13:28 PM
__author__ = "Michael Cuoco"

import logging, re

logger = logging.getLogger(__name__)  # configure logging
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


def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq: http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    array = array.flatten()  # all values are treated equally, arrays must be 1d
    if np.amin(array) < 0:
        array -= np.amin(array)  # values cannot be negative
    array += 0.0000001  # values cannot be 0
    array = np.sort(array)  # values must be sorted
    index = np.arange(1, array.shape[0] + 1)  # index per array element
    n = array.shape[0]  # number of array elements
    return (np.sum((2 * index - n - 1) * array)) / (
        n * np.sum(array)
    )  # Gini coefficient


TAGS = [
    "alignment_score",
    "L1_alignment_score",
    "mate_alignment_score",
    "alignment_score_normed",
    "L1_alignment_score_normed",
    "mate_alignment_score_normed",
    "L1_mapping_quality",
    "L1_reference_start",
    "L1_reference_end",
    "L1_Acount",
    "num_supp_alignments",
]


def features(p: dict) -> dict:
    """
    Extract features from a window of reads
    All reads must be read2 or a contig
    :param p: dict of region, must contain "reads" key
    """

    for n in ["Chromosome", "Start", "End", "reads"]:
        assert n in p, f"Region must contain {n}"

    # initialize lists for features
    l = dict()
    for k in TAGS + ["3end", "5end", "mapq", "starts"]:
        l[k] = []

    # initialize features dict
    f = {
        "Chromosome": p["Chromosome"],
        "Start": p["Start"],
        "End": p["End"],
        "n_reads": 0,
        "n_contigs": 0,
        "n_fwd": 0,
        "n_rev": 0,
        "n_duplicates": 0,
        "n_proper_pairs": 0,
        "n_ref_reads": 0,
        "max_mapq": 0,
        "min_mapq": 0,
        "n_unique_5end": 0,
        "n_unique_3end": 0,
        "3end_gini": float(0),
        "5end_gini": float(0),
    }

    for tag in TAGS:
        for n in [0, 0.25, 0.5, 0.75, 1]:
            f[tag + "_q" + str(n)] = float(0)
        f[tag + "_mean"] = float(0)

    # collect features from the reads in the window
    for i, r in enumerate(p["reads"]):
        if type(r) != Read:
            raise Exception("Reads must be of type Read")
        if r.is_read2:
            raise Exception("Reads must not be read2")

        if r.is_duplicate:
            f["n_duplicates"] += 1
            continue

        l["3end"].append(r.three_end)
        l["5end"].append(r.five_end)
        l["mapq"].append(r.mapping_quality)
        f["n_proper_pairs"] += r.is_proper_pair
        f["n_ref_reads"] += r.is_ref_read
        f["n_contigs"] += not r.is_read1

        if r.is_reverse:
            f["n_rev"] += 1
        else:
            f["n_fwd"] += 1

        for tag in TAGS:
            if "_normed" in tag:
                if getattr(r, tag.replace("_normed", "")):
                    if tag in [
                        "L1_alignment_score_normed",
                        "mate_alignment_score_normed",
                    ]:  # adjust alignments scores for read length
                        if not r.is_read1:
                            l[tag].append(
                                getattr(r, tag.replace("_normed", ""))
                                / getattr(r, "read_length")
                            )
                        else:
                            l[tag].append(
                                getattr(r, tag.replace("_normed", ""))
                                / getattr(r, "mate_read_length")
                            )
                    elif tag in ["alignment_score_normed"]:
                        l[tag].append(
                            getattr(r, tag.replace("_normed", ""))
                            / getattr(r, "read_length")
                        )
            elif getattr(r, tag):
                l[tag].append(getattr(r, tag))

    # compute mean and quantiles for these features
    for tag in TAGS:
        if len(l[tag]) > 0:
            quantiles = np.quantile(l[tag], [0, 0.25, 0.5, 0.75, 1])
            for n, q in zip([0, 0.25, 0.5, 0.75, 1], quantiles):
                f[tag + "_q" + str(n)] = float(q)
            f[tag + "_mean"] = np.mean(l[tag])

    f["n_reads"] = f["n_fwd"] + f["n_rev"]

    # if reads are empty, return empty features
    if f["n_reads"] == 0:
        return f

    f["3end_gini"] = gini(np.array(l["3end"], dtype=np.float64))
    f["5end_gini"] = gini(np.array(l["5end"], dtype=np.float64))
    f["max_mapq"] = max(l["mapq"])
    f["min_mapq"] = min(l["mapq"])
    f["n_unique_5end"] = len(set(l["5end"]))
    f["n_unique_3end"] = len(set(l["3end"]))

    return f
