#!/usr/bin/env python
# Created on: Nov 4, 2024 at 4:42:57 PM
__author__ = "Michael Cuoco"

from itertools import product
import numpy as np
import pandas as pd
import pyranges as pr
from time import time
import matplotlib.pyplot as plt
import seaborn as sns
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

HUE_ORDER = ["KNRGL", "OTHER", "KRGL"]


def _collate_labels(row):
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
    elif hasattr(row, "graffite") and row.graffite:
        return "KNRGL"
    elif hasattr(row, "xtea") and row.xtea:
        return "KNRGL"
    elif hasattr(row, "KNRGL") and row.KNRGL:
        return "KNRGL"
    else:
        return "OTHER"


def compute_distance(peaks, reference, name) -> pd.DataFrame:
    """
    Compute distance between peaks and germline annotations
    :param peaks: DataFrame of peaks
    :param reference: DataFrame of reference annotations
    :param name: Name of annotation
    """

    start = time()
    peaks.set_index("locus", inplace=True, drop=False)
    peak_cns = peaks[["Chromosome", "Start", "End", "locus"]].drop_duplicates()
    peak_cns = pr.PyRanges(peak_cns)
    reference.set_index("locus", inplace=True, drop=False)
    ref_cns = reference[["Chromosome", "Start", "End", "locus"]].drop_duplicates()
    ref_cns = pr.PyRanges(ref_cns)

    peak_distance = (
        peak_cns.nearest(ref_cns, overlap=False)
        .df.rename(columns={"Distance": name})
        .set_index("locus", drop=False)
    )
    peaks.loc[peak_distance.index, name] = peak_distance[name]

    return peaks.reset_index(drop=True)


def peak_stats(
    df: pd.DataFrame, by_label: bool = False, plot: bool = False
) -> pd.DataFrame:
    """
    Show distributions of some peak statistics
    :param df: DataFrame of peaks
    :param by_label: Whether to group by label
    :param plot: Whether to plot
    """

    # check columns
    stats = {
        "n_reads": True,
        "width": True,
        "germline_distance": True,
        "peak_distance": True,
    }
    for s in stats.copy():
        if s not in df.columns:
            stats.pop(s)
    print(f"Computing stats for {', '.join(stats.keys())}")

    # compute stats
    kwargs = {"kind": "ecdf", "col_wrap": 2, "log_scale": True}
    if by_label:
        out = df.groupby("label")[list(stats.keys())].describe().T
        df = df.melt(id_vars="label", value_vars=stats.keys(), var_name="statistic")
        kwargs["hue"] = "label"
        kwargs["hue_order"] = HUE_ORDER
    else:
        out = df[stats.keys()].describe()
        df = df.melt(value_vars=stats.keys(), var_name="statistic")

    # make plot
    if plot:
        sns.displot(df, x="value", col="statistic", **kwargs)

    return out


def label(peaks: pd.DataFrame, anno: pd.DataFrame, name: str) -> pd.DataFrame:
    cns = peaks[["Chromosome", "Start", "End", "locus"]].drop_duplicates()
    cns, anno = pr.PyRanges(cns), pr.PyRanges(anno)
    cns = cns.count_overlaps(anno, overlap_col=name)
    cns = cns.df.set_index("locus", drop=False)
    peaks = peaks.set_index("locus", drop=False)
    peaks.loc[cns.index, name] = cns[name]

    if peaks[name].isna().sum() > 0:
        raise ValueError(f"Failed to label {name}. Found NaNs")
    return peaks.reset_index(drop=True)


def label_peaks(peaks: pd.DataFrame, labels: dict[str]) -> pd.DataFrame:
    """
    Label peaks based on overlap with annotation files
    :param peaks: DataFrame of peaks
    :param labels: Dictionary of labels and their file paths
    """
    for c in ["Chromosome", "Start", "End", "locus"]:
        assert c in peaks.columns, f"Column '{c}' not in peaks DataFrame"

    # ensure locus is index of peaks
    peaks.set_index("locus", inplace=True, drop=False)

    print("Labeling peaks...")
    cns = peaks[["Chromosome", "Start", "End", "locus"]].drop_duplicates()
    germline = []
    for l, f in labels.items():
        anno = pr.read_bed(f).df
        germline.append(anno)
        cns = label(cns, anno, l)

        count = (cns[l] > 0).sum()
        print(f"Labelled {count} {l}")

        peaks.loc[cns["locus"], l] = cns[l]

    peaks["label"] = peaks.apply(_collate_labels, axis=1)

    return peaks.reset_index(drop=True)


def cells_per_peak(df) -> pd.Series:

    arrays = product([1, 5, 10, 50], ["n_cells", "n_regions", "n_donors"])
    index = pd.MultiIndex.from_tuples(arrays, names=["reads", "metric"])

    res = []
    for r in [1, 5, 10, 50]:
        mask = df["n_reads"].ge(r)
        for f in ["cell_id", "region", "donor_id"]:
            n = df.loc[mask, f].nunique()
            if f == "region":
                n = str(n)
            res.append(n)

    return pd.Series(res, index=index)


def _cluster_stats(df: pd.DataFrame) -> pd.Series:
    """
    Compute cluster statistics on a DataFrame of peaks
    """

    n_read_quantiles = np.quantile(df["n_reads"], [0, 0.25, 0.5, 0.75, 1])
    n_read_quantiles = {
        f"n_reads_q{int(v)}": q for v, q in zip([0, 25, 50, 75, 100], n_read_quantiles)
    }
    n_knrgl, n_other, n_krgl = (
        df["label"].value_counts().reindex(HUE_ORDER, fill_value=0)
    )
    n_reads = df["n_reads"].sum()
    if n_knrgl > 0:
        label = "KNRGL"
    elif n_krgl > 0:
        label = "KRGL"
    else:
        label = "OTHER"
    chrom, start, end = df["Chromosome"].iloc[0], df["Start"].min(), df["End"].max()
    return pd.Series(
        {
            "Chromosome": chrom,
            "Start": start,
            "End": end,
            "locus": f"{chrom}:{start}-{end}",
            "width": end - start,
            "n_peaks": df.shape[0],
            "n_KNRGL": n_knrgl,
            "n_KRGL": n_krgl,
            "n_OTHER": n_other,
            "n_reads": n_reads,
            "label": label,
            **n_read_quantiles,
        }
    )


def cluster_peaks(peaks, slack: int = 0) -> pd.DataFrame:

    for c in ["Chromosome", "Start", "End", "locus", "width", "n_reads", "label"]:
        assert c in peaks.columns, f"Column '{c}' not in peaks DataFrame"

    if "Cluster" in peaks.columns:
        print("Removing existing cluster column")
        peaks.drop(columns="Cluster", inplace=True)

    peaks.set_index("locus", inplace=True, drop=False)
    cns = peaks[["Chromosome", "Start", "End", "locus"]].drop_duplicates()
    cns = pr.PyRanges(cns).cluster(slack=slack).df
    cns.set_index("locus", inplace=True, drop=False)
    peaks.loc[cns.index, "Cluster"] = cns["Cluster"]
    clusters = peaks.groupby("Cluster").apply(_cluster_stats)
    print(f"Clustered {peaks.shape[0]:,} peaks into {clusters.shape[0]:,} clusters")
    print(clusters["label"].value_counts())

    return peaks, clusters
