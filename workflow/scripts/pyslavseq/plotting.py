#!/usr/bin/env python
# Created on: Apr 2, 2024 at 4:29:45 PM
__author__ = "Michael Cuoco"

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

HUE_ORDER = ["KNRGL", "OTHER", "KRGL"]


def peaks_per_cell(df, by_label=False):
    """
    Compute number of peaks per cell and plot
    :param df: DataFrame of peaks
    :param by_label: Group by label
    """
    for c in ["cell_id", "donor_id"]:
        assert c in df.columns, f"Column '{c}' not in DataFrame"

    g, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12))
    if not by_label:
        ppc = df.groupby(["cell_id", "donor_id"]).size().reset_index(name="n_peaks")
        total, mean, sd = (
            ppc["n_peaks"].sum(),
            ppc["n_peaks"].mean(),
            ppc["n_peaks"].std(),
        )
        sns.histplot(ppc, x="n_peaks", ax=ax1, bins=100)
        ax1.annotate(
            f"Total: {total:,}, Mean: {mean:.2f}+/-{sd:.2f}",
            xy=(0.65, 0.95),
            xycoords="axes fraction",
        )
        ax1.set_ylabel("Number of cells")
        sns.boxplot(ppc, x="donor_id", y="n_peaks", ax=ax2)
    else:
        ppc = (
            df.groupby(["cell_id", "donor_id", "label"])
            .size()
            .reset_index(name="n_peaks")
        )
        total = ppc.groupby("label")["n_peaks"].sum()
        mean = ppc.groupby("label")["n_peaks"].mean()
        sd = ppc.groupby("label")["n_peaks"].std()
        sns.histplot(
            ppc, x="n_peaks", hue="label", ax=ax1, bins=100, hue_order=HUE_ORDER
        )
        ax1.set_ylabel("Number of cells")
        for i, l in enumerate(HUE_ORDER):
            adj = (i + 1) * 0.05
            ax1.annotate(
                f"{l}: Total: {total[l]:,}, Mean: {mean[l]:.2f}+/-{sd[l]:.2f}",
                xy=(0.5, 1 - adj),
                xycoords="axes fraction",
            )
        sns.boxplot(
            ppc, x="donor_id", y="n_peaks", hue="label", ax=ax2, hue_order=HUE_ORDER
        )
        ax2.get_legend().remove()
    # rotate xlabels
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, horizontalalignment="right")


def datashader_plot(
    df: pd.DataFrame, x: str, y: str, log_scale=(False, False), **kwargs
):
    """
    Plot an x vs. y variable for all peaks, split by their labels in df using datashader

    TODO:
     - [ ]
    """

    df = df.copy()  # Avoid modifying the original DataFrame

    for c in [x, y, "label"]:
        assert c in df.columns, f"Column '{c}' not in DataFrame"

    if log_scale[0]:
        df[x] = np.log10(df[x])
    if log_scale[1]:
        df[y] = np.log10(df[y])

    nlabels = df["label"].nunique()
    fig, axs = plt.subplots(
        1, nlabels, figsize=(nlabels * 7, 6), sharey=True, sharex=True
    )
    if nlabels == 1:
        axs = [axs]
    fig.subplots_adjust(wspace=0)

    for ax, (l, d) in zip(axs, df.groupby("label")):
        artist = dsshow(
            d,
            ds.Point(x, y),
            ds.count(),
            norm="log",
            aspect="auto",
            x_range=(0, d[x].max()),
            y_range=(0, d[y].max()),
            ax=ax,
            **kwargs,
        )
        plt.colorbar(artist, ax=ax)

        ax.set_title(f"label = {l} (n = {d.shape[0]})")
        ax.set_xlabel(x)
        if log_scale[0]:
            ax.xaxis.set_major_formatter(
                FuncFormatter(lambda val, pos: f"$10^{{{int(val)}}}$")
            )
        if log_scale[1]:
            ax.yaxis.set_major_formatter(
                FuncFormatter(lambda val, pos: f"$10^{{{int(val)}}}$")
            )

    axs[0].set_ylabel(y)
