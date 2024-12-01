#!/usr/bin/env python
# Created on: Nov 7, 2024 at 1:02:36 PM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)

import pyBigWig
import pyranges as pr
import pandas as pd
import numpy as np

if __name__ == "__main__":

    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.DEBUG,
    )
    logger = logging.getLogger(__name__)

    peaks = pr.read_bed(snakemake.input.peaks).df  # type: ignore
    out = []
    # Group peaks by chromosome for efficient processing
    with pyBigWig.open(snakemake.input.bw) as bw:  # type: ignore
        for chrom, chrom_peaks in peaks.groupby("Chromosome", observed=True):
            if chrom not in bw.chroms():
                logger.warning(
                    f"Skipping chromosome {chrom}. No reads on this chromosome"
                )
                continue

            # Get all values for this chromosome at once
            logger.debug(f"Reading values for chromosome {chrom}")
            chrom_size = bw.chroms()[chrom]
            chrom_values = bw.values(chrom, 0, chrom_size, numpy=True)

            # Process each peak in this chromosome
            for peak in chrom_peaks.itertuples():
                # Extract values for this peak
                peak_values = chrom_values[peak.Start : peak.End]
                peak_values = peak_values[~np.isnan(peak_values)]

                if peak_values.shape[0] == 0:
                    logger.warning(
                        f"Skipping peak {chrom}:{peak.Start}-{peak.End}. No reads in this peak"
                    )
                    continue

                # Find maximum position and value
                max_idx = np.argmax(peak_values)
                max_val = peak_values[max_idx]
                abs_pos = peak.Start + max_idx

                logger.debug(
                    f"For peak {chrom}:{peak.Start}-{peak.End}, max {max_val} reads at {abs_pos + 1}"
                )

                out.append(
                    {
                        "Chromosome": chrom,
                        "Start": abs_pos + 1,  # Convert to 1-based coordinates
                        "End": abs_pos + 2,
                        "Score": max_val,
                    }
                )

    out = pd.DataFrame(out)
    pr.PyRanges(out).extend(250).df.to_csv(snakemake.output.bed, sep="\t", index=False, header=False)  # type: ignore
