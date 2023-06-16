#!/usr/bin/env python
# Created on: Jun 15, 2023 at 7:19:19 PM
__author__ = "Michael Cuoco"

import sys
from collections import defaultdict

import pandas as pd
import pyranges as pr
from pysam import VariantFile


def read_xtea_vcf(filename: str, rm_orphan=True):
    """Reads a XTEA VCF file and returns a pandas DataFrame."""
    out = defaultdict(list)
    for record in VariantFile(filename).fetch():

        # remove orphan calls if rm_orphan is True
        # orphan calls have higher FP rate
        if rm_orphan and "orphan" in record.info["SUBTYPE"]:
            continue

        out["Chromosome"].append(record.chrom)
        out["Start"].append(record.start)
        out["End"].append(record.stop)
        for k in record.info.keys():
            out[k].append(record.info[k])

    return pd.DataFrame(out)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    knrgl = read_xtea_vcf(snakemake.input[0])
    knrgl["repStart"] = knrgl["Start"]
    knrgl["repEnd"] = knrgl["End"]

    # save to BED
    pr.PyRanges(knrgl).to_bed(snakemake.output.knrgl)

    # save to BED with 1kb extension of 3end
    knrgl_1kb_3nd = knrgl.copy()
    knrgl_1kb_3nd["Start"] = knrgl_1kb_3nd["Start"].apply(
        lambda x: x["Start"] - 1000 if x["Strand"] == "+" else x["Start"], axis=1
    )
    knrgl_1kb_3nd["End"] = knrgl_1kb_3nd["End"].apply(
        lambda x: x["End"] - 1000 if x["Strand"] == "-" else x["End"], axis=1
    )
    knrgl_1kb_3nd["Start"] = knrgl_1kb_3nd["Start"].apply(lambda x: 0 if x < 0 else x)
    pr.PyRanges(knrgl_1kb_3nd).to_bed(snakemake.output.knrgl_1kb_3nd)

    # save to BED with 20kb extensions of both ends
    knrgl_20kb = knrgl.copy()
    knrgl_20kb["Start"] -= 20000
    knrgl_20kb["End"] += 20000
    knrgl_20kb["Start"] = knrgl_20kb["Start"].apply(lambda x: 0 if x < 0 else x)
    pr.PyRanges(knrgl_20kb).to_bed(snakemake.output.knrgl_20kb)

    sys.stderr.close()
