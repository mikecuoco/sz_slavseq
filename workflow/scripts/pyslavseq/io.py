#!/usr/bin/env python
# Created on: Nov 4, 2024 at 1:11:51 PM
__author__ = "Michael Cuoco"

# configure logging
import logging

logger = logging.getLogger(__name__)

from tqdm import tqdm
from collections import defaultdict
from time import time
from pathlib import Path
import re
import pandas as pd
import pyranges as pr
from functools import cache


def write_tabix(df: pd.DataFrame, path: str):
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


def read_tabix(fn: str) -> pd.DataFrame:

    from pyarrow import csv

    meta_fn = "/iblm/netapp/data3/mcuoco/sz_slavseq/config/all_samples.tsv"
    meta = pd.read_csv(
        meta_fn, sep="\t", dtype={"sample_id": str, "tissue_id": str, "donor_id": str}
    )
    meta = meta[["sample_id", "tissue_id", "region"]].set_index("sample_id")

    print("Reading peaks from " + fn)
    start = time()
    df = csv.read_csv(fn, parse_options=csv.ParseOptions(delimiter="\t")).to_pandas()
    df.columns = [c.lstrip("#") for c in df.columns]
    df["tissue_id"] = df["cell_id"].map(meta["tissue_id"])
    df["region"] = df["cell_id"].map(meta["region"])
    if "locus" not in df.columns:
        df["locus"] = (
            df["Chromosome"]
            + ":"
            + df["Start"].astype(str)
            + "-"
            + df["End"].astype(str)
        )
    print(
        f"Read {len(df):,} peaks from {df.cell_id.nunique()} cells in {time() - start:.2f} seconds"
    )

    return df


def read_megane_mei(f: str) -> pd.DataFrame:
    """
    Read MEGaNe genotyped MEI bed file output
    :param f: Path to MEGaNe MEI bed file
    """

    def parse_mei(string):
        match = re.search(r"ref_pos=(\d+),chimeric=(\d+),hybrid=(\d+)", string)
        if match:
            return {
                "ref_pos": int(match.group(1)),
                "chimeric": int(match.group(2)),
                "hybrid": int(match.group(3)),
            }

    def parse_confidence(string):
        return string.split(":")[1] if ":" in string else None

    def parse_unique(string):
        match = re.search(r"unique:(\w+),50bp_or_longer:(\w+),orig_conf:(\w+)", string)
        if match:
            return {
                "unique_status": match.group(1),
                "longer_50bp": match.group(2),
                "orig_conf": match.group(3),
            }

    def parse_subfamily(string):
        match = re.search(r"status=(\w+),MEI=(\w+)", string)
        if match:
            return {"subfamily_status": match.group(1), "subfamily_MEI": match.group(2)}

    def parse_genotype_quality(string):
        match = re.search(r"(\d);(\w+);(\w)", string)
        if match:
            return {"genotype": int(match.group(1)), "quality": match.group(2)}

    def parse_genotype_TSD(string):
        match = re.search(r"(\w+);(\d+);(\w+)", string)
        if match:
            return {"TSD": match.group(3)}

    # Define the column names and converters for each column to parse
    column_names = [
        "Chromosome",
        "Start",
        "End",
        "MEI",
        "MEI_left",
        "MEI_right",
        "confidence",
        "unique",
        "subfamily_pred",
        "transduction",
        "genotype_quality",
        "genotype_TSD",
        "bp_left",
        "bp_right",
        "ID",
    ]

    converters = {
        "MEI_left": parse_mei,
        "MEI_right": parse_mei,
        "confidence": parse_confidence,
        "unique": parse_unique,
        "subfamily_pred": parse_subfamily,
        "genotype_quality": parse_genotype_quality,
        "genotype_TSD": parse_genotype_TSD,
    }

    df = pd.read_csv(f, sep="\t", names=column_names, converters=converters)

    for col in [
        "MEI_left",
        "MEI_right",
        "unique",
        "subfamily_pred",
        "genotype_quality",
        "genotype_TSD",
    ]:
        expanded_df = df[col].apply(pd.Series)
        if "MEI_" in col:
            expanded_df = expanded_df.add_prefix(f"{col}_")
        df = df.drop(columns=[col]).join(expanded_df)

    return df


def read_megane_mea(f: str) -> pd.DataFrame:
    """
    Read MEGaNe genotyped MEA bed file output
    :param f: Path to MEGaNe MEA bed file
    """

    def parse_repeat(value):
        """Parse the 'repeat' column to split by ':'."""
        return {"subfamily_ME": value.split(":")[0], "ME": value.split(":")[1]}

    def parse_sequences(value):
        """Parse the 'sequences' column to split by ';'."""
        return value.split(";")

    def parse_genotype_filter(value):
        """Parse the 'pass_status' column to split by ';' and handle multiple status flags."""
        return {"genotype": int(value.split(";")[0]), "filter": value.split(";")[1]}

    def extract_tsd_length(value):
        """Extract numeric TSD length from 'TSD_len' column."""
        return float(value.split("=")[1])

    # Define converters for specific columns, excluding the 'sequences' column
    converters = {
        "repeat": parse_repeat,
        "genotype_filter": parse_genotype_filter,
        "TSD_len": extract_tsd_length,
    }

    # Load data with specified converters and exclude the 'sequences' column
    df = pd.read_csv(
        f,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3, 5, 6],
        names=[
            "Chromosome",
            "Start",
            "End",
            "repeat",
            "TSD_len",
            "genotype_filter",
        ],
        converters=converters,
    )
    df[["genotype", "filter"]] = df["genotype_filter"].apply(pd.Series)
    df.drop(columns=["genotype_filter"], inplace=True)
    df[["subfamily_ME", "ME"]] = df["repeat"].apply(pd.Series)
    df.drop(columns=["repeat"], inplace=True)

    return df


def read_rmsk_bed(file: str):
    "read rmsk bed file into dataframe"

    coord_conv = lambda x: int(x.rstrip(")").lstrip("("))

    return pd.read_csv(
        file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
        names=[
            "Chromosome",
            "Start",
            "End",
            "repName",
            "Score",
            "Strand",
            "milliDiv",
            "milliDel",
            "milliIns",
            "genoLeft",
            "repClassFamily",
            "repStart",
            "repEnd",
            "repLeft",
        ],
        converters={
            "repStart": coord_conv,
            "repLeft": coord_conv,
            "genoLeft": coord_conv,
        },
    )


@cache
def load_rmsk():
    """
    Load and process RepeatMasker annotations if not already cached
    """

    print("Loading RepeatMasker annotations")
    rmsk = read_rmsk_bed(
        "/iblm/netapp/data3/mcuoco/sz_slavseq/resources/chm13v2.0.XY/chm13v2.0.XY.fasta.all_rmsk.bed"
    )
    rmsk["locus"] = (
        rmsk["Chromosome"]
        + ":"
        + rmsk["Start"].astype(str)
        + "-"
        + rmsk["End"].astype(str)
    )
    my_rmsk = {}
    for i in ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]:
        my_rmsk[i] = rmsk.query("repName == @i and repEnd > 6000").reset_index()
        my_rmsk[i] = pr.PyRanges(my_rmsk[i]).three_end().extend({"5": 200, "3": 100}).df

    for i in ["(A)n", "(T)n"]:
        my_rmsk[i] = rmsk.query("repName == @i").reset_index(drop=True)

    return my_rmsk


@cache
def _load_megane_files():
    """
    Load and process MEGaNe annotations if not already cached
    """
    print("Loading MEGaNe annotations")
    METADATA = "/iblm/netapp/data3/mcuoco/sz_slavseq/config/all_donors.tsv"
    libd2donor = (
        pd.read_csv(METADATA, sep="\t")[["libd_id", "donor_id"]]
        .set_index("libd_id")
        .squeeze()
        .to_dict()
    )
    WGSDIR = (
        "/iblm/netapp/data3/mcuoco/sz_slavseq/resources/chm13v2.0.XY/wgs_calls/30x/"
    )
    megane_files = {}
    for l, d in libd2donor.items():
        megane_files[("MEI", d)] = Path(WGSDIR) / l / "MEI_final_gaussian_genotyped.bed"
        megane_files[("MEA", d)] = Path(WGSDIR) / l / "absent_MEs_genotyped.bed"

    return megane_files


def load_donor_genotype(donor: str) -> dict:
    """
    Load genotyped MEI and MEA calls for a donor
    """

    # Get cached rmsk data or load if not cached
    my_rmsk = load_rmsk()

    # Get cached MEGaNe files or load if not cached
    megane_files = _load_megane_files()

    # get nonreference MEI calls
    mei = (
        read_megane_mei(megane_files[("MEI", donor)])
        .query("confidence == 'high;PASS' and subfamily_MEI == 'L1HS'")
        .reset_index()
    )
    mei = pr.PyRanges(mei).extend(10).df
    mei_het = mei.query("genotype == 1").reset_index()
    mei_hom = mei.query("genotype == 2").reset_index()

    # get reference MEA calls
    mea = (
        read_megane_mea(megane_files[("MEA", donor)])
        .query('filter == "PASS" and subfamily_ME == "L1HS"')
        .reset_index()
    )
    l1hs = my_rmsk["L1HS"]
    mea_het = mea.query("genotype == 1").reset_index()
    l1hs = (
        pr.PyRanges(l1hs).count_overlaps(pr.PyRanges(mea_het), overlap_col="mea_het").df
    )
    mea_hom = mea.query("genotype == 2").reset_index()
    l1hs = (
        pr.PyRanges(l1hs).count_overlaps(pr.PyRanges(mea_hom), overlap_col="mea_hom").df
    )

    return {
        "MEI_het": mei_het,
        "MEI_hom": mei_hom,
        "L1HS_hom": l1hs.query("mea_hom == 0 and mea_het == 0"),
        "L1HS_het": l1hs.query("mea_hom == 0 and mea_het > 0"),
    }
