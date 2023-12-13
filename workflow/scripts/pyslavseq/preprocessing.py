#!/usr/bin/env python
# Created on: Aug 8, 2023 at 1:53:30 PM
__author__ = "Michael Cuoco"

# source: Genome In A Bottle (GIAB) genome stratifications
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/

# configure logging
import logging

logger = logging.getLogger(__name__)

import pandas as pd

STRATIFICATIONS = {
    "chm13": {
        "nonunique": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/CHM13@all/Mappability/CHM13_nonunique_l100_m2_e1.bed.gz",
        "difficult": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/validation/CHM13@all/Union_CHM13_alldifficultregions.bed.gz",
        "benchmark": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/CHM13v2.0_HG002-T2TQ100-V1.0_stvar.benchmark.bed",
    },
    "grch38": {
        "nonunique": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/Mappability/GRCh38_nonunique_l100_m2_e1.bed.gz",
        "difficult": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz",
        "benchmark": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed",
    },
}


def get_stratification(genome: str, name: str) -> pd.DataFrame:
    """
    Get specified genome stratification for grch38 or chm13 from GIAB FTP
    """

    if genome not in STRATIFICATIONS:
        raise Exception(
            f"Invalid genome {genome}, must be one of {list(STRATIFICATIONS.keys())}"
        )
    if name not in STRATIFICATIONS[genome]:
        raise Exception(
            f"Invalid stratification {name}, must be one of {list(STRATIFICATIONS[genome].keys())}"
        )

    return pd.read_csv(
        STRATIFICATIONS[genome][name],
        sep="\t",
        header=None,
        skiprows=1,
        names=["Chromosome", "Start", "End"],
    )
