#!/usr/bin/env python
# Created on: Apr 2, 2024 at 4:29:45 PM
__author__ = "Michael Cuoco"

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import datashader as ds
import datashader.transfer_functions as tf


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


def germline_dist(df: pd.DataFrame, y: str, hue=None):
    """
    Plot the germline distance vs. a given y variable for all peaks in df.
    Use log scale for x-axis. Use datashaer for visualization.
    """

    assert "germline_dist" in df.columns, "Column 'germline_dist' not in DataFrame"
    assert y in df.columns, f"Column '{y}' not in DataFrame"

    if not hue:
        df["germline_dist_log"] = np.log10(df["germline_dist"])
        nlabels = df["label"].nunique()
        g, axs = plt.subplots(
            1, nlabels, figsize=(5 * nlabels, 5), sharey=True, sharex=True
        )
        g.subplots_adjust(wspace=0)
        cvs = ds.Canvas(
            plot_width=400,
            plot_height=400,
            x_range=(0, df["germline_dist_log"].max()),
            y_range=(df[y].min(), df[y].max()),
        )

        for ax, (l, d) in zip(axs, df.groupby("label")):
            agg = cvs.points(d, "germline_dist_log", y)
            img = tf.shade(agg, how="log").to_pil()
            ax.imshow(
                img,
                extent=[0, df["germline_dist_log"].max(), df[y].min(), df[y].max()],
                aspect="auto",
            )
            ax.set_title(f"label = {l}")
            ax.set_xlabel("Distance to nearest germline (bp)")
            ax.set_xlim(0, df["germline_dist_log"].max())
            ax.xaxis.set_major_formatter(
                FuncFormatter(lambda val, pos: f"$10^{{{int(val)}}}$")
            )

        axs[0].set_ylabel(y)

    else:
        assert hue in df.columns, f"Column '{hue}' not in DataFrame"
        g = sns.relplot(
            data=df,
            x="germline_dist",
            y=y,
            hue=hue,
            kind="scatter",
            col="label",
            alpha=0.5,
            s=3,
        )
        g.set(xscale="log")

    return g
