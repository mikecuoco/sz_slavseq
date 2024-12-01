#!/usr/bin/env python
# Created on: Nov 4, 2024 at 5:57:55 PM
__author__ = "Michael Cuoco"

import numpy as np
import pandas as pd
import pyranges as pr

# compute cluster stats
def get_cluster_stats(df, peaks):
    cells = peaks.loc[df.name]
    n_read_quantiles = np.quantile(cells["n_reads"], [0, 0.25, 0.5, 0.75, 1])
    n_read_quantiles = {
        f"n_reads_q{int(v)}": q for v, q in zip([0, 25, 50, 75, 100], n_read_quantiles)
    }
    n_knrgl, n_other, n_krgl = (
        df["label"].value_counts().reindex(HUE_ORDER, fill_value=0)
    )
    if n_krgl > 0:
        label = "KRGL"
    elif n_knrgl > 0:
        label = "KNRGL"
    else:
        label = "OTHER"
    return pd.Series(
        {
            "Chromosome": df["Chromosome"].iloc[0],
            "Start": df["Start"].min(),
            "End": df["End"].max(),
            "width": df["End"].max() - df["Start"].min(),
            "germline_distance": df["germline_distance"].max(),
            "peak_distance": df["peak_distance"].max(),
            "n_peaks": df.shape[0],
            "n_KNRGL": n_knrgl,
            "n_KRGL": n_krgl,
            "n_OTHER": n_other,
            "label": label,
            **n_read_quantiles,
        }
    )


def cluster_peaks(peaks, slack: int = 0) -> pd.DataFrame:

    for c in [
        "Chromosome",
        "Start",
        "End",
        "locus",
        "width",
        "n_reads",
        "label",
        "germline_distance",
        "peak_distance",
    ]:
        assert c in peaks.columns, f"Column '{c}' not in peaks DataFrame"

    peaks = pr.PyRanges(peaks).cluster(slack=slack).df
    clusters = peaks.groupby("Cluster").apply(get_cluster_stats, peaks=peaks)

    return clusters
