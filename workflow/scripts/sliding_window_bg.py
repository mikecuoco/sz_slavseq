#!/usr/bin/env python
# Created on: Oct 21, 2024 at 7:25:38 AM
__author__ = "Michael Cuoco"

import logging, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import time
from functools import partial

# data science
from scipy.stats import poisson, false_discovery_control
import numpy as np
import pandas as pd

# bio
import pyBigWig

logger = logging.getLogger(__name__)


def process_chromosome(
    bw: pyBigWig.pyBigWig, chrom: str, wsize: int, bgsizes: list, minreads: int
) -> pd.DataFrame:
    """Process a single chromosome to calculate signal enrichment over background.

    Args:
            bw: Open pyBigWig file handle
            chrom: Chromosome name
            wsize: Window size
            bgsizes: List of background window sizes
            minreads: Minimum number of reads required
    """
    logger.debug(f"Reading values for {chrom}")
    start = time()
    sizes = [wsize] + bgsizes
    half_window = wsize // 2
    chrom_size = bw.chroms(chrom)

    # Get values and create position array
    values = pd.Series(bw.values(chrom, 0, chrom_size, numpy=True)).fillna(0)
    pos = np.arange(chrom_size)
    logger.debug(f"Read {len(values):,} values for {chrom}")

    # Calculate rolling means for all window sizes
    logger.debug(f"Calculating means for {chrom}")
    means = np.column_stack([values.rolling(s, center=True).mean() for s in sizes])
    means = np.nan_to_num(means, nan=0, copy=False)

    # Calculate rolling maxes for the window size for minread filtering
    wmax = values.rolling(wsize, center=True).max().fillna(0).values
    del values

    # Apply minreads and nonzero masks
    minreads_mask = wmax >= minreads
    del wmax
    nonzero_mask = ~(means[:, 1:] == 0).all(axis=1)
    mask = minreads_mask | nonzero_mask
    means = means[mask]
    pos = pos[mask]
    bgmax = means[:, 1:].max(axis=1)

    return pos, means[:, 0], bgmax


def sliding_window_bg(
    bw: str, wsize: int, bgsizes: list, minreads: int, n_threads: int = 1
) -> pd.DataFrame:
    """Calculate signal enrichment windows and background windows using threaded processing.

    Args:
            bw: Path to BigWig file
            wsize: Window size
            bgsizes: List of background window sizes
            minreads: Minimum reads required
            n_threads: Number of threads to use
    """
    half_window = wsize // 2
    bgsizes = list(bgsizes)

    results = {
        "Chromosome": [],
        "pos": [],
        "mean": [],
        "bgmax": [],
    }

    with pyBigWig.open(bw) as bw:
        chroms = [c for c in bw.chroms() if c != "chrM"]

        # Process chromosomes in parallel
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            process_func = partial(
                process_chromosome, bw, wsize=wsize, bgsizes=bgsizes, minreads=minreads
            )
            futures = {executor.submit(process_func, chrom): chrom for chrom in chroms}
            for future in as_completed(futures):
                chrom = futures[future]
                pos, means, bgmax = future.result()

                # append results
                results["Chromosome"].extend([chrom] * len(pos))
                results["pos"].append(pos)
                results["mean"].append(means)
                results["bgmax"].append(bgmax)
                logger.debug(f"Finished {chrom} in {time() - start:.2f} seconds")

    # Combine and process results
    results["pos"] = np.concatenate(results["pos"], dtype=np.int32)
    results["mean"] = np.concatenate(results["mean"], dtype=np.float32)
    results["bgmax"] = np.concatenate(results["bgmax"], dtype=np.float32)

    return results


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

    with pyBigWig.open(snakemake.input.bw) as bw:
        header = [(c, s) for c, s in bw.chroms().items() if c != "chrM"]

    windows = sliding_window_bg(
        snakemake.input.bw,
        wsize=snakemake.params.wsize,
        bgsizes=snakemake.params.bgsizes,
        minreads=snakemake.params.minreads,
        n_threads=snakemake.threads,
    )

    # Calculate derived values
    # logger.debug(f"Calculating statistics for {len(windows['pos']):,} loci...")
    # windows["bg_max"] = np.max(windows["mean"], axis=1)
    # windows["fc"] = windows["mean"][:, 0] / windows["bg_max"]
    # windows["pval"] = np.log10(poisson._sf(windows["mean"], windows["bg_max"]), dtype=np.float32)*-1
    # windows["qval"] = np.log10(false_discovery_control(windows["pval"], method="bh"), dtype=np.float32)*-1

    # cleanup memory
    # del windows["mean"], windows["bg_max"]

    # write to bigwig
    outfiles = {
        "wmeans": snakemake.output.wmeans,
        "bgmax": snakemake.output.bgmax,
    }
    for i, (signal, fname) in enumerate(outfiles.items()):
        with pyBigWig.open(fname, "w") as bw:
            bw.addHeader(header)

            # Write chromosome by chromosome
            for chrom, _ in header:
                chrom_mask = np.array(
                    [x == chrom for x in windows["Chromosome"]], dtype=bool
                )
                positions = windows["pos"][chrom_mask]
                values = windows["mean"][chrom_mask][:, i]

                logger.debug(
                    f"Writing {len(values):,} {signal} mean values for {chrom}"
                )
                bw.addEntries(str(chrom), positions, values=values, span=1)
