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
        if rm_orphan and ("orphan" in record.info["SUBTYPE"]):
            continue

        out["Chromosome"].append(record.chrom)
        out["Start"].append(record.start)
        out["End"].append(record.stop)
        for k in record.info.keys():
            out[k].append(record.info[k])

    return pd.DataFrame(out).rename(columns={"STRAND": "Strand"})


def extend_3end(x, extend=1000):
    if x["Strand"] == "+":
        x["End"] += extend
    elif x["Strand"] == "-":
        x["Start"] -= extend
    else:
        x["Start"] -= extend / 2
        x["End"] += extend / 2

    if x["Start"] < 0:
        x["Start"] = 0

    return x


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    knrgl = read_xtea_vcf(snakemake.input[0])
    knrgl["repStart"] = knrgl["Start"]
    knrgl["repEnd"] = knrgl["End"]

    # save to BED
    knrgl.to_csv(snakemake.output.knrgl, sep="\t", index=False)

    # save to BED with 1kb extension of 3end
    knrgl_1kb_3end = knrgl.copy()
    knrgl_1kb_3end = knrgl_1kb_3end.apply(extend_3end, axis=1)
    pr.PyRanges(knrgl_1kb_3end).sort().df.to_csv(
        snakemake.output.knrgl_1kb_3end, sep="\t", index=False
    )

    # save to BED with 20kb extensions of both ends
    knrgl_20kb = knrgl.copy()
    knrgl_20kb["Start"] -= 2e4
    knrgl_20kb["End"] += 2e4
    knrgl_20kb.loc[knrgl_20kb["Start"] < 0, "Start"] = 0
    pr.PyRanges(knrgl_20kb).sort().df.to_csv(
        snakemake.output.knrgl_20kb, sep="\t", index=False
    )

    sys.stderr.close()
