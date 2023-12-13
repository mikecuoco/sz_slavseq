#!/usr/bin/env python
# Created on: Nov 21, 2023 at 12:13:28 PM
__author__ = "Michael Cuoco"

import logging, re

logger = logging.getLogger(__name__)  # configure logging
from collections import namedtuple
from pysam import AlignedSegment
import numpy as np


def is_ref_read(read: AlignedSegment) -> bool:
    "return True if read is ref, False if non-ref based on mate tag"

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
        end = cigar[-1] if read.is_reverse else cigar[0]
        clipped = int(end[0]) if end[1] == "H" or end[1] == "S" else 0
        return clipped < 30


# use named tuple to convert AlignedSegment
# gives a speedup over using pysam.AlignedSegment (about 1.5x)
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
        "is_secondary",
        "is_supplementary",
        "is_mapped",
        "mapping_quality",
        "num_supp_alignments",
        "cigarstring",
        "alignment_score",
        "L1_alignment_score",
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
        read.reference_name,
        read.reference_start,
        read.reference_end,
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
        is_ref_read(read),
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


def basic_stats(p: dict):
    """
    Calculate basic stats about a region
    :param p: dict of region, must contain "reads" key
    """

    for n in ["Chromosome", "Start", "End", "reads"]:
        assert n in p, f"Region must contain {n}"

    # iterate over reads to calculate stats
    for mate in ["r1", "r2"]:
        for k in ["n_ref_reads", "max_mapq", "n_duplicates", "n_fwd", "n_rev"]:
            p[f"{mate}_{k}"] = 0
        p[f"{mate}_starts"] = []

    for r in p["reads"]:
        mate = "r1" if r.is_read1 else "r2"
        if r.is_duplicate:
            p[f"{mate}_n_duplicates"] += 1
            continue

        if r.is_reverse:
            p[f"{mate}_n_rev"] += 1
        else:
            p[f"{mate}_n_fwd"] += 1

        p[f"{mate}_max_mapq"] = max(p[f"{mate}_max_mapq"], r.mapping_quality)
        p[f"{mate}_n_ref_reads"] += r.is_ref_read

        if r.is_forward:
            p[f"{mate}_starts"].append(r.reference_start)
        else:
            p[f"{mate}_starts"].append(r.reference_end)

    for mate in ["r1", "r2"]:
        p[f"{mate}_n_unique_starts"] = len(set(p[f"{mate}_starts"]))
        p[f"{mate}_n_reads"] = p[f"{mate}_n_fwd"] + p[f"{mate}_n_rev"]
        del p[f"{mate}_starts"]
    p["n_reads"] = p["r1_n_reads"] + p["r2_n_reads"]

    # add strand info
    if (p["r1_n_fwd"] == 0) and (p["r1_n_rev"] > 0):
        p["Strand"] = "-"
    elif (p["r1_n_rev"] == 0) and (p["r1_n_fwd"] > 0):
        p["Strand"] = "+"
    else:
        p["Strand"] = "."

    # remove reads
    del p["reads"]

    return p


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


def features(p: dict) -> dict:
    """
    Extract features from a window of reads
    All reads must be read1
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
        "n_fwd": 0,
        "n_rev": 0,
        "n_duplicates": 0,
        "n_proper_pairs": 0,
        "n_ref_reads": 0,
        "3end_gini": float(0),
        "5end_gini": float(0),
        "max_mapq": 0,
        "orientation_bias": float(0),
        "frac_proper_pairs": float(0),
        "frac_duplicates": float(0),
        "n_unique_starts": 0,
    }

    for tag in TAGS:
        for n in [0, 0.25, 0.5, 0.75, 1]:
            f[tag + "_q" + str(n)] = float(0)
        f[tag + "_mean"] = float(0)

    # collect features from the reads in the window
    for i, r in enumerate(p["reads"]):
        if type(r) != Read:
            raise Exception("Reads must be of type Read")
        if not r.is_read1:
            raise Exception("Reads must all be read1")

        if i == 0:
            start = r.reference_start

        if r.is_duplicate:
            f["n_duplicates"] += 1
            continue

        l["starts"].append(r.reference_start)
        if not r.is_reverse:
            l["3end"].append(r.reference_end - start)
            l["5end"].append(r.reference_start - start)
        else:
            l["3end"].append(r.reference_start - start)
            l["5end"].append(r.reference_end - start)

        l["mapq"].append(r.mapping_quality)
        f["n_proper_pairs"] += r.is_proper_pair
        f["n_ref_reads"] += r.is_ref_read

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
                f[tag + "_q" + str(n)] = q
            f[tag + "_mean"] = np.mean(l[tag])

    f["n_reads"] = f["n_fwd"] + f["n_rev"]

    # if reads are empty, return empty features
    if f["n_reads"] == 0:
        return f

    f["3end_gini"] = gini(np.array(l["3end"], dtype=np.float64))
    f["5end_gini"] = gini(np.array(l["5end"], dtype=np.float64))
    f["max_mapq"] = max(l["mapq"])
    f["n_unique_starts"] = len(set(l["starts"]))

    return f
