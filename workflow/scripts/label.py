#!/usr/bin/env python
# Created on: Nov 17, 2023 at 6:45:57 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pyranges as pr
from pysam import VariantFile
import pyBigWig


def read_megane_LINE1(filename: str, libd_id: str):
    """
    Read MEGAnE LINE1 calls into a pandas dataframe
    Filter for LINE/L1 insertions with FILTER=PASS
    """

    var_df = defaultdict(list)
    for rec in VariantFile(filename).fetch():
        if rec.info["SVTYPE"][0] == "LINE/L1":
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
    var_df = (
        pd.DataFrame(var_df)
        .explode("individual")
        .query("individual == @libd_id")
        .reset_index(drop=True)
    )

    logger.info(f"Loaded {len(var_df)} L1HS insertions for {libd_id} from {filename}")

    return var_df


def read_xtea_LINE1(filename: str, libd_id: str):
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

    var_df = (
        pd.DataFrame(var_df)
        .query("~SUBTYPE.str.contains('transduction')")
        .explode("individual")
        .query("individual == @libd_id")
        .reset_index(drop=True)
    )

    logger.info(
        f"Loaded {len(var_df)} non-transduction insertions for {libd_id} from {filename}"
    )

    return var_df


def label_nearest(
    x: pd.DataFrame,
    y: pd.DataFrame,
    name: str,
    slack: int = 0,
) -> pd.DataFrame:
    """
    Label windows in x that are nearest to windows in y

    Parameters
    ----------
    x : pd.DataFrame, windows to label
    y : pd.DataFrame, windows to overlap with
    name : str, name of column to add to df
    slack : int, slack for pyranges join
    """

    if "cell_id" in x.columns:
        assert (
            x["cell_id"].nunique() == 1
        ), "Cannot label windows with multiple cell ids"

    y[f"{name}_id"] = y.index.values
    y = y[[f"{name}_id", "Chromosome", "Start", "End"]]

    nearest = (
        pr.PyRanges(y)
        .nearest(pr.PyRanges(x))
        .df.query("Distance <= @slack")
        .drop(columns=["Distance", "Start", "End"])
        .rename(columns={"Start_b": "Start", "End_b": "End"})
    )

    merged = x.merge(nearest, on=["Chromosome", "Start", "End"], how="left")
    x[name] = merged[f"{name}_id"].notnull()
    x[f"{name}_id"] = merged[f"{name}_id"]

    # check if any labels are duplicated
    assert (
        x[x[name]][f"{name}_id"].duplicated().sum() == 0
    ), f"Duplicate labels for {name}"

    return x


def label_overlap(x: pd.DataFrame, y: pd.DataFrame, name: str):
    x = pr.PyRanges(x).count_overlaps(pr.PyRanges(y), overlap_col=name).df
    x[name] = x[name].astype(bool)
    return x


def add_final_features(df: pd.DataFrame, en_pos_bw, en_neg_bw):
    """
    Add features to dataframe
    """

    # add en_motif scores
    df["en_pos_score"] = df.apply(
        lambda x: np.max(en_pos_bw.values(x["Chromosome"], x["Start"], x["End"])),
        axis=1,
    )
    df["en_neg_score"] = df.apply(
        lambda x: np.max(en_neg_bw.values(x["Chromosome"], x["Start"], x["End"])),
        axis=1,
    )

    # add final features
    df["width"] = df["End"] - df["Start"]
    df["rpm"] = df["n_reads"] * (df["size_factor"] / 1e6)
    df["frac_contigs"] = df["n_contigs"] / df["n_reads"]
    df["orientation_bias"] = np.maximum(df["n_fwd"], df["n_rev"]) / df["n_reads"]
    df["frac_proper_pairs"] = df["n_proper_pairs"] / df["n_reads"]
    df["frac_duplicates"] = df["n_duplicates"] / (df["n_reads"] + df["n_duplicates"])
    df["frac_unique_3end"] = df["n_unique_3end"] / df["n_reads"]
    df["frac_unique_5end"] = df["n_unique_5end"] / df["n_reads"]
    df["frac_mean_supp_alignments"] = df["num_supp_alignments_mean"] / df["n_reads"]

    return df


def label_ref_peaks(df: pd.DataFrame):
    """
    Using n_ref_reads, find peaks across all cells from a single donor
    """

    df = pr.PyRanges(df).cluster().df
    labelled = df.groupby("Cluster").apply(lambda d: d["n_ref_reads"].sum() > 0)
    labelled.name = "ref"
    return df.set_index("Cluster").join(labelled).reset_index()


if __name__ == "__main__":
    from pyslavseq.preprocessing import read_rmsk_bed

    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )

    logger = logging.getLogger(__name__)

    # read in annotations
    rmsk = read_rmsk_bed(snakemake.input.rmsk)  # type: ignore
    LINE1 = [
        "L1HS_3end",
        "L1PA2_3end",
        "L1PA3_3end",
        "L1PA4_3end",
        "L1PA5_3end",
        "L1PA6_3end",
    ]
    rmsk = rmsk[(rmsk["repName"].isin(LINE1)) & (rmsk["repEnd"] > 860)]

    indv_to_libd = (
        pd.read_csv(snakemake.config["donors"], sep="\t")
        .set_index("donor_id")["libd_id"]
        .to_dict()
    )  # type: ignore
    libd_id = indv_to_libd[int(snakemake.wildcards.donor)]  # type: ignore

    annotations = {
        "xtea": read_xtea_LINE1(snakemake.input.xtea_vcf, libd_id),
        "megane_perc": read_megane_LINE1(
            snakemake.input.megane_percentile_vcf, libd_id
        ),
        "megane_gaus": read_megane_LINE1(snakemake.input.megane_gaussian_vcf, libd_id),
        "primer_sites": pr.read_bed(snakemake.input.primer_sites).df,
    }

    for l1 in LINE1:
        annotations[l1] = rmsk.query("repName == @l1").reset_index(drop=True)

    # bulk slavseq
    bulk = pq.read_table(snakemake.input.bulk).to_pandas().reset_index(drop=True)
    bulk["REF"] = bulk["n_ref_reads"] > 0
    annotations["bulk_NR"] = bulk.query("REF == False").reset_index(drop=True)
    annotations["bulk_R"] = bulk.query("REF == True").reset_index(drop=True)

    # open bigwigs
    en_pos_bw = pyBigWig.open(snakemake.input.en_motif_pos, "r")  # type: ignore
    en_neg_bw = pyBigWig.open(snakemake.input.en_motif_neg, "r")  # type: ignore

    logger.info(f"Reading regions for donor {snakemake.wildcards.donor}..")  # type: ignore
    data = []
    for f in snakemake.input.cells:  # type: ignore
        logger.info(f"Reading {f}")
        cell = pq.read_table(f).to_pandas().reset_index(drop=True)

        # add cell id
        cell["cell_id"] = Path(f).name.rstrip(".pqt")
        cell["donor_id"] = snakemake.wildcards.donor

        # annotate
        for name, anno in annotations.items():
            logger.info(f"Labelling {name}")
            cell = label(cell, anno, name, 200)

        # add final features
        cell = add_final_features(cell, en_pos_bw, en_neg_bw)

        data.append(cell)

    data = pd.concat(data)
    logger.info("Labelling ref peaks")
    data = label_ref_peaks(data)
    data.to_parquet(snakemake.output[0])  # type: ignore

    en_pos_bw.close()
    en_neg_bw.close()
