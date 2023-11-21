#!/usr/bin/env python
# Created on: Jun 15, 2023 at 7:30:33 PM
__author__ = "Michael Cuoco"

import sys
from myutils.rmsk import read_rmsk
import pyranges as pr


def adjust_coordinates(x, adjust=5900):
    "Adjust coordinates of L1s to be >= 5900bp from the 3'-end"

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
    "Parse RepeatMasker output file and return a DataFrame of L1s with intact 3' ends."

    rmsk = read_rmsk(filename)

    # keep L1s in these 6 families with intact 3' ends
    rep_names = [
        "L1HS_3end",
        "L1PA2",
        "L1PA3",
        "L1PA4",
        "L1PA5",
        "L1PA6",
        "L1PA7",
        "L1PA8",
    ]
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


def fix_negative_ends(df):
    # TODO: take min() of coordinates that are larger than chromosome length
    df["Start"] = df.Start.apply(lambda x: max(0, x))
    df["End"] = df.End.apply(lambda x: max(0, x))
    return df


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    rmsk = read_rmsk_L1_3ends(snakemake.input[0], min_repend=5000)
    rmsk["repStart"] = rmsk["Start"]
    rmsk["repEnd"] = rmsk["End"]

    # save to BED
    pr.PyRanges(rmsk).to_bed(snakemake.output.rmsk)

    # save to BED with 1kb extension of 3end
    rmsk_1kb_3end = rmsk.copy()
    rmsk_1kb_3end = pr.PyRanges(rmsk_1kb_3end).extend({"3": 1000}).sort().df
    rmsk_1kb_3end = fix_negative_ends(rmsk_1kb_3end)
    pr.PyRanges(rmsk_1kb_3end).to_bed(snakemake.output.rmsk_1kb_3end)

    sys.stderr.close()
