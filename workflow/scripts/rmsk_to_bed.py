#!/usr/bin/env python
# Created on: Jun 15, 2023 at 7:30:33 PM
__author__ = "Michael Cuoco"

import sys
from myutils.rmsk import read_rmsk
import pyranges as pr


def adjust_coordinates(x, adjust=5900):
    if x["strand"] == "-":
        if x["repLeft"] < adjust:
            x["genoEnd"] -= adjust - x["repLeft"]
    else:
        if x["repStart"] < adjust:
            x["genoStart"] += adjust - x["repStart"]

    if x["genoEnd"] <= x["genoStart"]:
        x["genoEnd"] = x["genoStart"] + 1

    assert (
        x["genoEnd"] > x["genoStart"]
    ), f"genoEnd <= genoStart for {x['genoEnd']} <= {x['genoStart']}"
    assert x["genoStart"] >= 0, f"genoStart < 0 for {x['genoStart']} < 0"
    return x


def read_rmsk_L1_3ends(filename: str, min_repend: int = 5900):
    "Parse RepeatMasker output file and return a DataFrame of L1s intact 3' ends."

    rmsk = read_rmsk(filename)

    # keep L1s in these 6 families with intact 3' ends
    rep_names = ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]
    rmsk = rmsk.loc[
        (rmsk["repName"].isin(rep_names)) & (rmsk["repEnd"] > min_repend), :
    ]

    rmsk = rmsk.apply(adjust_coordinates, axis=1)
    rmsk = rmsk.rename(
        columns={
            "genoName": "Chromosome",
            "genoStart": "Start",
            "genoEnd": "End",
            "strand": "Strand",
            "repName": "Name",
        }
    )
    rmsk = rmsk.loc[
        (rmsk.Start >= 0) & (rmsk.End >= 0),
        ["Chromosome", "Start", "End", "Strand", "Name"],
    ]

    return rmsk


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    rmsk = read_rmsk_L1_3ends(snakemake.input[0])
    rmsk["repStart"] = rmsk["Start"]
    rmsk["repEnd"] = rmsk["End"]

    # save to BED
    pr.PyRanges(rmsk).to_bed(snakemake.output.rmsk)

    # save to BED with 1kb extension of 3end
    rmsk_1kb_3end = rmsk.copy()
    pr.PyRanges(rmsk_1kb_3end).extend({"3": 1000}).sort().to_bed(
        snakemake.output.rmsk_1kb_3end
    )

    # save to BED with 20kb extensions of both ends
    rmsk_20kb = rmsk.copy()
    pr.PyRanges(rmsk_20kb).extend({"3": 2e4, "5": 2e4}).sort().to_bed(
        snakemake.output.rmsk_20kb
    )

    sys.stderr.close()
