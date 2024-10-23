#!/usr/bin/env python
# Created on: Oct 21, 2024 at 7:25:38 AM
__author__ = "Michael Cuoco"

import logging
from time import time
import numpy as np
import pandas as pd
import pyBigWig
import bionumpy as bnp
from bionumpy.arithmetics import merge_intervals
from bionumpy.datatypes import Interval
from pandas.api.indexers import BaseIndexer

logger = logging.getLogger(__name__)


class GenomicIndexer(BaseIndexer):
    def __init__(self, genomic_positions, window_size):
        # Ensure genomic positions are provided as a pandas Series
        self.genomic_positions = pd.Series(genomic_positions)
        self.window_size = window_size

    def get_window_bounds(self, num_values, min_periods, center, closed, step):
        # Create arrays to store the start and end of each window
        start = pd.Series(0, index=self.genomic_positions.index)
        end = pd.Series(0, index=self.genomic_positions.index)

        # Iterate over genomic positions and define the windows
        for i in range(num_values):
            # Define the start and end of the window based on genomic distances
            current_pos = self.genomic_positions.iloc[i]
            start_idx = max(
                0, self.genomic_positions.searchsorted(current_pos - self.window_size)
            )
            end_idx = self.genomic_positions.searchsorted(
                current_pos + self.window_size, side="right"
            )

            start.iloc[i] = start_idx
            end.iloc[i] = end_idx

        return start.to_numpy(), end.to_numpy()


def sliding_window(
    bw: str, size: int = 200, step: int = 1, minreads: int = 5
) -> pd.DataFrame:
    """
    Find peaks in a bigwig file using a sliding window approach
    :param bw: path to bigwig file
    :param size: size of the window
    :param step: step size for the window
    :param minreads: minimum number of reads in a window to be considered a peak
    :return: pandas DataFrame of peaks
    """

    # TODO: add background enrichment test?

    with pyBigWig.open(bw) as b:
        out = []
        for chrom, chrom_size in b.chroms().items():
            start = time()
            logger.info(f"Finding peaks for {chrom}")
            values = b.values(
                chrom, 0, chrom_size, numpy=True
            )  # get values for chromosome from bigwig
            values = pd.Series(values).dropna()  # remove zeros
            logger.debug("calculating rolling max")
            indexer = GenomicIndexer(values.index, window_size=size)
            maxes = values.rolling(indexer, step=step, min_periods=1).max()
            logger.debug(f"keep windows with >= {minreads}...")
            starts = maxes.index[maxes >= minreads]
            if len(starts) == 0:
                logger.warning(f"no peaks found for {chrom}")
                continue
            intervals = [(chrom, start, start + 1) for start in starts]
            logger.info(
                f"Found {len(intervals)} windows with >= {minreads} reads in {time() - start:.2f} seconds"
            )
            peaks = merge_intervals(
                Interval.from_entry_tuples(intervals), distance=size
            )
            logger.info(f"Merged into {len(peaks)} peaks")
            out.append(peaks.topandas())

    return pd.concat(out)


if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )
    logger = logging.getLogger(__name__)
    logger.info(f"Numpy version: {np.__version__}")
    logger.info(f"Pandas version: {pd.__version__}")
    logger.info(f"PyBigWig version: {pyBigWig.__version__}")
    logger.info(f"bionumpy version: {bnp.__version__}")

    df = sliding_window(
        snakemake.input.bw,
        size=snakemake.params.size,
        step=snakemake.params.step,
        minreads=snakemake.params.minreads,
    )  # type: ignore

    # write
    df.to_csv(snakemake.output.bed, sep="\t", index=False, header=False)  # type: ignore
