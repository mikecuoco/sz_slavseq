#!/usr/bin/env python
# Created on: Apr 2, 2024 at 4:29:45 PM
__author__ = "Michael Cuoco"

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import datashader as ds
from datashader.mpl_ext import dsshow
from mpl_toolkits.axes_grid1 import ImageGrid


def joint_ecdfplot(df: pd.DataFrame, x: str, y: str, **kwargs):
    """
    Plot the number of cells and donors per peak in a joint plot.
    :param df: DataFrame with columns "n_cells" and a specified y variable.
    :param kwargs: Additional keyword arguments to pass to seaborn functions.
    """

    for c in [x, y]:
        assert c in df.columns, f"Column '{c}' not in DataFrame"
    if "hue" in kwargs:
        assert kwargs["hue"] in df.columns, f"Column '{kwargs['hue']}' not in DataFrame"

    g = sns.JointGrid(data=df, x=x, y=y)
    sns.ecdfplot(data=df, x=x, ax=g.ax_marg_x, **kwargs)
    sns.ecdfplot(data=df, y=y, ax=g.ax_marg_y, **kwargs)

    g.plot_joint(sns.scatterplot, data=df, alpha=0.3, s=3, **kwargs)
    if g.ax_marg_x.get_legend():
        g.ax_marg_y.get_legend().remove()
        g.ax_joint.get_legend().remove()

    return g


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
