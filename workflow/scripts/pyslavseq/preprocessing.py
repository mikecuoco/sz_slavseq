#!/usr/bin/env python
# Created on: Aug 8, 2023 at 1:53:30 PM
__author__ = "Michael Cuoco"

# configure logging
import logging

logger = logging.getLogger(__name__)

import pandas as pd
import pyranges as pr


def collate_labels(row):
    if hasattr(row, "n_ref_reads") and row.n_ref_reads > 0:
        return "KRGL"
    elif hasattr(row, "primer_sites") and row.primer_sites:
        return "KRGL"
    elif hasattr(row, "rmsk") and row.rmsk:
        return "KRGL"
    elif hasattr(row, "ref") and row.ref:
        return "KRGL"
    elif hasattr(row, "bulk_ref") and row.bulk_ref:
        return "KRGL"
    elif hasattr(row, "l1hs") and row.l1hs:
        return "KRGL"
    # elif hasattr(row, "l1pa2") and row.l1pa2:
    #     return "KRGL"
    # elif hasattr(row, "l1pa3") and row.l1pa3:
    #     return "KRGL"
    # elif hasattr(row, "l1pa4") and row.l1pa4:
    #     return "KRGL"
    # elif hasattr(row, "l1pa5") and row.l1pa5:
    #     return "KRGL"
    # elif hasattr(row, "l1pa6") and row.l1pa6:
    #     return "KRGL"
    elif hasattr(row, "megane") and row.megane:
        return "KNRGL"
    # elif hasattr(row, "megane_breakpoints") and row.megane_breakpoints:
    #     return "KNRGL"
    elif hasattr(row, "graffite") and row.graffite:
        return "KNRGL"
    elif hasattr(row, "xtea") and row.xtea:
        return "KNRGL"
    elif hasattr(row, "KNRGL") and row.KNRGL:
        return "KNRGL"
    else:
        return "OTHER"


# source: Genome In A Bottle (GIAB) genome stratifications
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/

STRATIFICATIONS = {
    "chm13": {
        "nonunique": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/CHM13@all/Mappability/CHM13_nonunique_l100_m2_e1.bed.gz",
        "difficult": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/validation/CHM13@all/Union_CHM13_alldifficultregions.bed.gz",
        "benchmark": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107/CHM13v2.0_HG002-T2TQ100-V1.0_stvar.benchmark.bed",
        "telomere": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_telomere.bed",
        "censat": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.0.bed",
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


def df2tabix(df: pd.DataFrame, path: str):
    """
    Write a pandas DataFrame to a tabix-indexed file
    """
    import pysam

    assert path.endswith(".bed.gz"), "Path must end with .bed"
    if not df.columns[0].startswith("#"):
        print("Renaming " + df.columns[0] + " to #" + df.columns[0])
        df.columns = ["#" + c if i == 0 else c for i, c in enumerate(df.columns)]

    for c, d in zip(["#Chromosome", "Start", "End"], df.columns[0:2]):
        assert c == d, f"Column {c} not found in DataFrame but is required"

    # write to file
    df.to_csv(path.rstrip(".gz"), sep="\t", index=False)

    # create tabix index
    pysam.tabix_index(path.rstrip(".gz"), preset="bed", force=True)

    return path
