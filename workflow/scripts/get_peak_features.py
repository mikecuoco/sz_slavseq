import logging


logger = logging.getLogger(__name__)

from pathlib import Path
import pysam
import pyranges as pr
import numpy as np
import pandas as pd
from make_bw_peaks import GenomicIndexer
from time import time


def gini(array):
    """
    Calculate the Gini coefficient of a numpy array
    from https://github.com/oliviaguest/gini
    """
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


TAGS = {
    "three_end_clippedA": "CA",
    "three_end_clipped_length": "CL",
    "three_end_clippedA_normed": "CA",
    "three_end_clipped_length_normed": "CL",
    "L1_mapping_quality": "LQ",
    "L1_alignment_score": "L1",
    "L1_alignment_score_normed": "L1",
    "L1_Acount": "LA",
    "L1_Acount_normed": "LA",
    "mate_alignment_score": "MS",
    "mate_alignment_score_normed": "MS",
}

FEATURES = [
    "n_reads",
    "n_contigs",
    "n_fwd",
    "n_rev",
    "n_duplicates",
    "n_proper_pairs",
    "n_ref_reads",
    "n_unique_5end",
    "n_unique_3end",
    "n_unique_clipped_3end",
    "clipped_3end_gini",
    "3end_gini",
    "5end_gini",
    *[
        f"{t}_mean" for t in list(TAGS.keys()) + ["mapq"]
    ],  # this might fail because of the TAGS variable
    *[
        f"{t}_q{n}"
        for t in list(TAGS.keys()) + ["mapq"]
        for n in [0, 0.25, 0.5, 0.75, 1]
    ],
]

# TODO: compute breakpoints?
def features(p: dict) -> dict:
    """
    Extract features from a peak of reads
    All reads must be read1 or a contig
    :param p: dict of region, must contain "reads" key
    """
    for key in ["Chromosome", "Start", "End", "reads"]:
        assert key in p, f"Peak must contain {key}"

    # initialize lists for features
    l = {k: [] for k in list(TAGS.keys()) + ["3end", "5end", "mapq", "clipped_3end"]}
    f = {
        "Chromosome": p["Chromosome"],
        "Start": p["Start"],
        "End": p["End"],
        **{k: 0 for k in FEATURES},
    }

    # collect features from the reads in the peak
    logger.info(
        f"Extracting features from {len(p['reads'])} for {p['Chromosome']}:{p['Start']}-{p['End']}"
    )
    start = time()
    for r in p["reads"]:
        # check if read is valid
        if not isinstance(r, pysam.AlignedSegment):
            logger.error(f"Read is not of type pysam.AlignedSegment: {r}")
            raise ValueError("Reads must be of type pysam.AlignedSegment")
        if r.is_read2:
            logger.error(f"Read is read2: {r}")
            raise ValueError("Reads must not be read2")

        # if duplicate, add to count and skip
        if r.is_duplicate:
            f["n_duplicates"] += 1
            continue

        l["3end"].append(r.reference_start if r.is_reverse else r.reference_end)
        l["5end"].append(r.reference_end if r.is_reverse else r.reference_start)
        l["mapq"].append(r.mapping_quality)
        f["n_proper_pairs"] += r.is_proper_pair
        f["n_ref_reads"] += r.get_tag("RR")
        f["n_contigs"] += not r.is_read1
        f["n_reads"] += 1
        f["n_rev" if r.is_reverse else "n_fwd"] += 1

        if r.get_tag(TAGS["three_end_clipped_length"]) > 0:
            l["clipped_3end"].append(
                r.reference_start if r.is_reverse else r.reference_end
            )

        # add read length based features
        read_length = r.infer_read_length()
        mate_read_length = r.get_tag("ML") if r.is_read1 else read_length

        for name, tag in TAGS.items():
            if "_normed" in tag:
                if r.has_tag(tag.replace("_normed", "")):
                    norm_tag = r.get_tag(tag.replace("_normed", ""))
                    if tag in [
                        "L1_alignment_score_normed",
                        "mate_alignment_score_normed",
                    ]:
                        l[name].append(norm_tag / mate_read_length)
                    elif tag == "alignment_score_normed":
                        l[name].append(norm_tag / read_length)
            elif r.has_tag(tag):
                l[name].append(r.get_tag(tag))

    assert (
        len(l["L1_Acount"]) == f["n_reads"]
    ), "3end length does not match number of reads"

    logger.info(f"Finished in {time() - start:.2f} seconds")

    # compute mean and quantiles for these features
    for tag in list(TAGS.keys()) + ["mapq"]:
        if l[tag]:
            quantiles = np.quantile(l[tag], [0, 0.25, 0.5, 0.75, 1])
            for n, q in zip([0, 0.25, 0.5, 0.75, 1], quantiles):
                f[f"{tag}_q{n}"] = float(q)
            f[f"{tag}_mean"] = np.mean(l[tag])

    # check if features are in correct order
    if list(f.keys()) != ["Chromosome", "Start", "End", *FEATURES]:
        logger.error(f"Features are not in correct order: {list(f.keys())}")
        raise ValueError("Features are not in correct order")

    # if empty, return empty features
    if f["n_reads"] == 0:
        logger.info(f"No reads in peak {p['Chromosome']}:{p['Start']}-{p['End']}")
        return f

    f["3end_gini"] = gini(np.array(l["3end"], dtype=np.float64))
    f["5end_gini"] = gini(np.array(l["5end"], dtype=np.float64))
    f["clipped_3end_gini"] = (
        gini(np.array(l["clipped_3end"], dtype=np.float64)) if l["clipped_3end"] else 0
    )
    f["n_unique_5end"], f["n_unique_3end"] = len(set(l["5end"])), len(set(l["3end"]))
    f["n_clipped_3end"], f["n_unique_clipped_3end"] = len(l["clipped_3end"]), len(
        set(l["clipped_3end"])
    )

    return f


if __name__ == "__main__":

    from pyslavseq.preprocessing import df2tabix

    logger = logging.getLogger(__name__)

    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )

    read_filter = (
        lambda r: (not r.is_read2)
        and r.is_mapped
        and (not r.is_supplementary)
        and (not r.is_secondary)
    )

    # read the peaks from the bed file
    logger.info(f"Reading peaks from {snakemake.input.bed}")
    peaks = pr.read_bed(snakemake.input.bed).df  # type: ignore

    # open the bam file, extract reads in the peaks, and extract features
    logger.info(
        f"Extracting features for {len(peaks)} peaks from reads of {snakemake.input.bam}"
    )
    with pysam.AlignmentFile(snakemake.input.bam, "rb") as bam:  # type: ignore
        res = []
        for p in peaks.itertuples():
            # get the reads in the peak from the bam
            reads = bam.fetch(p.Chromosome, p.Start, p.End)
            reads = filter(read_filter, reads)

            # save the reads to a dict
            w = {
                "Chromosome": p.Chromosome,
                "Start": p.Start,
                "End": p.End,
                "reads": list(reads),
            }

            # get the features for the peak
            res.append(features(w))

    peaks = pd.DataFrame.from_records(res).query("n_reads > 0").reset_index(drop=True)

    logger.info("Adding additional features")
    peaks["width"] = peaks["End"] - peaks["Start"]
    peaks["locus"] = (
        peaks["Chromosome"].astype(str)
        + ":"
        + peaks["Start"].astype(str)
        + "-"
        + peaks["End"].astype(str)
    )
    peaks["orientation_bias"] = (
        np.maximum(peaks["n_fwd"], peaks["n_rev"]) / peaks["n_reads"]
    )
    peaks["frac_proper_pairs"] = peaks["n_proper_pairs"] / peaks["n_reads"]
    peaks["frac_duplicates"] = peaks["n_duplicates"] / (
        peaks["n_reads"] + peaks["n_duplicates"]
    )
    peaks["frac_unique_3end"] = peaks["n_unique_3end"] / peaks["n_reads"]
    peaks["frac_unique_5end"] = peaks["n_unique_5end"] / peaks["n_reads"]
    peaks["frac_clipped_3end"] = peaks["n_clipped_3end"] / peaks["n_reads"]
    peaks["frac_unique_clipped_3end"] = (
        peaks["n_unique_clipped_3end"] / peaks["n_clipped_3end"]
    )

    # TODO compute features of neighboring peaks?

    cell_id = Path(snakemake.input.bam).name.rstrip(".tagged.sorted.bam")
    peaks["cell_id"] = cell_id
    donor_id = Path(snakemake.input.bam).parent.name
    peaks["donor_id"] = donor_id

    # write to file
    logger.info(f"Writing features to {snakemake.output.bed}")
    df2tabix(peaks, snakemake.output.bed)  # type: ignore
