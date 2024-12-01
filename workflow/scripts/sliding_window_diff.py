#!/usr/bin/env python
# Created on: Nov 30, 2024 at 11:17:52 PM
__author__ = "Michael Cuoco"

import logging
from time import time

# data science
import numpy as np
import pandas as pd

# bio
import pyranges as pr
import pyBigWig

logger = logging.getLogger(__name__)


def rolling_max(bw, chrom, chrom_size, wsize):
    """
    Calculate rolling maximum over a genomic window.

    Parameters:
            bw (pyBigWig.BigWig): BigWig file object
            chrom (str): Chromosome name
            chrom_size (int): Chromosome size
            wsize (int): Window size for rolling calculation

    Returns:
            pd.Series: Rolling maximum values
    """
    values = bw.values(chrom, 0, chrom_size)
    values = pd.Series(values).fillna(0)
    return values.rolling(wsize, center=True).max()


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log[0],
        filemode="w",
        level=logging.DEBUG,
    )
    logger = logging.getLogger(__name__)
    logger.debug(f"Numpy version: {np.__version__}")
    logger.debug(f"Pandas version: {pd.__version__}")
    logger.debug(f"PyBigWig version: {pyBigWig.__version__}")

    half_window = snakemake.params.wsize // 2

    hippo = pyBigWig.open(snakemake.input.hippo[0])
    dlpfc = pyBigWig.open(snakemake.input.dlpfc[0])
    for chrom, chrom_size in hippo.chroms().items():
        if chrom == "chrM":
            continue
        logger.debug(f"Processing {chrom}...")

        # Get values
        hippo_max = rolling_max(hippo, chrom, chrom_size, snakemake.params.wsize)
        dlpfc_max = rolling_max(dlpfc, chrom, chrom_size, snakemake.params.wsize)

        # Get background
        hippo_bg = rolling_max(hippo, chrom, chrom_size, snakemake.params.bgsize)
        dlpfc_bg = rolling_max(dlpfc, chrom, chrom_size, snakemake.params.bgsize)

        # Filter
        mask = hippo_bg.ge(snakemake.params.minreads) | dlpfc_bg.ge(
            snakemake.params.minreads
        )
        for v in [hippo_max, hippo_bg, dlpfc_max, dlpfc_bg]:
            v = v[mask]

        # Create DataFrame once with all data
        res = pd.DataFrame(
            {
                "Chromosome": chrom,
                "Start": hippo_max.index - half_window,
                "End": hippo_max.index + half_window,
                "hippo": hippo_max,
                "hippo_bg": hippo_bg,
                "dlpfc": dlpfc_max,
                "dlpfc_bg": dlpfc_bg,
            }
        )

        # Vectorized likelihood calculation
        start = time()
        # Add small epsilon to avoid log(0)
        eps = 1e-10
        res["loglr"] = (
            hippo_max * (np.log(hippo_max + eps) - np.log(dlpfc_max + eps))
            + dlpfc_max
            - hippo_max
        )
        logger.debug(
            f"Computed likelihood for {len(res):,} regions in {time() - start:.2f} seconds"
        )

        # Find exclusive windows
        mask = hippo_max.gt(0) & dlpfc_max.eq(0)
        hippo_only = res[mask].reset_index(drop=True)
        mask = hippo_max.eq(0) & dlpfc_max.gt(0)
        dlpfc_only = res[mask].reset_index(drop=True)
        logger.debug(
            f"Found {len(hippo_only)} hippo-only and {len(dlpfc_only)} dlpfc-only windows"
        )

        # TODO Write output

    # close files
    hippo.close()
    dlpfc.close()
