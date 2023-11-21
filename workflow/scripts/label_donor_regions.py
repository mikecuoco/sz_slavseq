#!/usr/bin/env python
# Created on: Nov 17, 2023 at 6:45:57 PM
__author__ = "Michael Cuoco"

import logging

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

logger = logging.getLogger(__name__)

import pyarrow.parquet as pq
import pyranges as pr
import pandas as pd
from collections import defaultdict
from pysam import VariantFile


def read_megane_LINE1(filename: str):
    """
    Read MEGAnE LINE1 calls into a pandas dataframe
    Filter for LINE/L1 insertions with FILTER=PASS
    """

    var_df = defaultdict(list)
    for rec in VariantFile(filename).fetch():
        if "LINE/L1" in rec.info["SVTYPE"]:
            var_df["Chromosome"].append(rec.chrom)
            var_df["Start"].append(rec.pos)
            var_df["End"].append(rec.stop)
            var_df["AC"].append(rec.info["AC"])
            var_df["MEI"].append(rec.info["MEI"][0])
            var_df["FILTER"].append(rec.filter.keys()[0])
            var_df["SVLEN"].append(rec.info["SVLEN"])
            var_df["MESTRAND"].append(rec.info["MESTRAND"][0])
            var_df["MEPRED"].append(rec.info["MEPRED"][0])
            indv = []
            for i in rec.samples.keys():
                if 1 in rec.samples[i]["GT"]:
                    indv.append(i)
            var_df["individual"].append(indv)

    var_df = pd.DataFrame(var_df)

    ins_per_indv = (
        var_df.query("FILTER == 'PASS'")
        .explode("individual")
        .groupby("individual")
        .count()["Chromosome"]
        .mean()
    )

    logger.info(
        f"Loaded {len(var_df)} variants from {filename} with {ins_per_indv:.2f} FILTER=PASS insertions per individual"
    )

    return var_df


def read_xtea_LINE1(filename: str):
    """
    Read XTEA LINE1 calls into a pandas dataframe
    Remove transduction insertions (SUBTYPE contains "transduction")
    """

    var_df = defaultdict(list)
    for rec in VariantFile(filename).fetch():
        var_df["Chromosome"].append(rec.chrom)
        var_df["Start"].append(rec.pos)
        var_df["End"].append(rec.stop)
        var_df["SUBTYPE"].append(rec.info["SUBTYPE"])
        var_df["SVLEN"].append(rec.info["SVLEN"])
        var_df["STRAND"].append(rec.info["STRAND"][0])
        indv = []
        for i in rec.samples.keys():
            if 1 in rec.samples[i]["GT"]:
                indv.append(i.rstrip(".md"))
        var_df["individual"].append(indv)
    var_df = pd.DataFrame(var_df)

    ins_per_indv = (
        var_df.query("~SUBTYPE.str.contains('transduction')")
        .explode("individual")
        .groupby("individual")
        .count()["Chromosome"]
        .mean()
    )

    logger.info(
        f"Loaded {len(var_df)} variants from {filename} with {ins_per_indv:.2f} non-transduction insertions per individual"
    )

    return var_df


def label(
    df: pd.DataFrame,
    other_df: pd.DataFrame,
    name: str,
    add_id: bool = False,
    slack: int = 0,
) -> pd.DataFrame:
    """
    Label windows in df with whether they overlap with other_df

    Parameters
    ----------
    df : pd.DataFrame, windows to label
    other_df : pd.DataFrame, windows to overlap with
    name : str, name of column to add to df
    add_id : bool, whether to add a column with the id of the overlapping window
    slack : int, slack for pyranges join
    """

    # add id column if requested
    if add_id:
        other_df.loc[:, f"{name}_id"] = other_df.index.values
        other_df = other_df[[f"{name}_id", "Chromosome", "Start", "End"]]
    else:
        other_df = other_df[["Chromosome", "Start", "End"]]

    # join with pyranges to find overlaps
    overlap = (
        pr.PyRanges(df[["Chromosome", "Start", "End"]], int64=True)
        .join(
            pr.PyRanges(other_df, int64=True),
            how="left",
            slack=slack,
        )
        .df
    )

    overlap[name] = overlap["Start_b"] != -1
    overlap = overlap.drop(columns=["Start_b", "End_b"]).drop_duplicates()

    # check for duplicated rows
    df = df.join(
        overlap.set_index(["Chromosome", "Start", "End"]),
        on=["Chromosome", "Start", "End"],
        how="left",
    )

    if (
        df.shape[0]
        != df[["Chromosome", "Start", "End", "cell_id"]].drop_duplicates().shape[0]
    ):
        logger.error(f"Some overlaps are duplicated for {name}")
        raise ValueError(f"Some overlaps are duplicated for {name}")

    return df


logger.info(f"Reading regions for donor {snakemake.wildcards.donor}..")  # type: ignore
data = pq.read_table(snakemake.input.regions).to_pandas()  # type: ignore
data = data.sort_values(["Chromosome", "Start", "End"])

# TODO: label difficult regions

logger.info(f"Read donor metadata from {snakemake.config['donors']}")  # type: ignore
meta = pd.read_csv(snakemake.config["donors"], sep="\t")  # type: ignore
indv_to_libd = meta.set_index("donor_id")["libd_id"].to_dict()
libd_id = indv_to_libd[int(snakemake.wildcards.donor)]  # type: ignore

logger.info(f"Labelling megane calls using {snakemake.input.megane_vcf}")  # type: ignore
megane = (
    read_megane_LINE1(snakemake.input.megane_vcf)  # type: ignore
    .query("Chromosome != 'chrY'")
    .query("FILTER == 'PASS'")
    .explode("individual")
    .query("individual == @libd_id")
)
data = label(data, megane, "megane")
data = label(data, megane, "megane_500bp", slack=500)
data = label(data, megane, "megane_1000bp", slack=1000)

logger.info(f"Labelling xtea calls from {snakemake.input.xtea_vcf}")  # type: ignore
xtea = (
    read_xtea_LINE1(snakemake.input.xtea_vcf)  # type: ignore
    .explode("individual")
    .query("~SUBTYPE.str.contains('transduction')")
    .query("individual == @libd_id")
)
data = label(data, xtea, "xtea")
data = label(data, xtea, "xtea_500bp", slack=500)
data = label(data, xtea, "xtea_1000bp", slack=1000)

logger.info(f"Labelling SLAVseq primer sites from {snakemake.input.primer_sites}")  # type: ignore
primer_sites = pr.read_bed(snakemake.input.primer_sites).df  # type: ignore
data = label(data, primer_sites, "primer_sites")

logger.info(f"Labelling donor bulk peaks from {snakemake.input.bulk_peaks}")  # type: ignore
bulk_peaks = pr.read_bed(snakemake.input.bulk_peaks).df  # type: ignore
data = label(data, bulk_peaks, "bulk_peaks")

data.to_parquet(snakemake.output[0])  # type: ignore
