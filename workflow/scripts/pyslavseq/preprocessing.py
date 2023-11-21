#!/usr/bin/env python
# Created on: Aug 8, 2023 at 1:53:30 PM
__author__ = "Michael Cuoco"

# configure logging
import logging

logger = logging.getLogger(__name__)

import pandas as pd


def get_hg38_blacklist() -> pd.DataFrame:
    """
    Get blacklist MHC, KIR, Tandem Repeat, SegDups, Gaps, and False duplicated regions from NCBI for hg38 genome builds
    """
    # read blacklist regions from NCBI
    region_urls = {
        "mhc": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_MHC.bed.gz",
        "kir": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_KIR.bed.gz",
        "trs": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeats_201to10000bp_slop5.bed.gz",
        "segdups": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/SegmentalDuplications/GRCh38_segdups.bed.gz",
        "gaps": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_gaps_slop15kb.bed.gz",
        "false_dup": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_false_duplications_correct_copy.bed.gz",
    }

    # add to dictionary
    regions = []
    for id, url in region_urls.items():
        regions.append(
            pd.read_csv(
                url,
                sep="\t",
                header=None,
                skiprows=1,
                names=["Chromosome", "Start", "End"],
            )
        )
        regions[-1]["Name"] = id

    return pd.concat(regions)
