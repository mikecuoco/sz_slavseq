#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import sys
from pathlib import Path
import pyarrow.parquet as pq
import pandas as pd
import pyranges as pr
from joblib import Parallel, delayed


def label_windows(df: pd.DataFrame, other_df: pd.DataFrame, name: str = "Name"):
    assert type(df) == pd.DataFrame, "df must be a pandas DataFrame"
    assert type(other_df) == pd.DataFrame, "other_df must be a pandas DataFrame"

    # join with pyranges to find overlaps
    overlap = (
        pr.PyRanges(df[["Chromosome", "Start", "End"]], int64=True)
        .join(
            pr.PyRanges(other_df[["Chromosome", "Start", "End"]], int64=True),
            how="left",
        )
        .df
    )

    overlap[name] = overlap["Start_b"] != -1
    overlap = overlap.drop(columns=["Start_b", "End_b"]).drop_duplicates()

    return df.join(
        overlap.set_index(["Chromosome", "Start", "End"]),
        on=["Chromosome", "Start", "End"],
        how="left",
    )


def collate_labels(x):

    # check if blacklist is in the series
    if x.isin(["blacklist"]).any() and x["blacklist"]:
        return "blacklist"

    for i in ["", "_1kb_3end"]:
        for l1 in ["xtea", "L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]:
            if x[l1 + i]:
                return l1 + i

    for l1 in ["xtea", "L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]:
        if x[l1 + "_20kb"]:
            return l1 + "_20kb"

    return "unknown"


def read_features(file: str, anno: dict) -> pd.DataFrame:

    df = pq.read_table(file).to_pandas()
    df = df.loc[df["n_ref_reads"] == 0, :]  # filter out reference windows
    # check if df is empty
    if df.empty:
        return pd.DataFrame()

    df["cell_id"] = Path(file).stem.rstrip("_windows")
    for l in anno.keys():
        df = label_windows(df, anno[l], l)
    df["label"] = df.apply(collate_labels, axis=1)

    return df


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    # read in the annotation files
    anno = {
        "xtea": pd.read_csv(snakemake.input.xtea, sep="\t"),
        "xtea_1kb_3end": pd.read_csv(snakemake.input.xtea_1kb_3end, sep="\t"),
        "xtea_20kb": pd.read_csv(snakemake.input.xtea_20kb, sep="\t"),
    }

    # add bulk peaks if available
    if snakemake.input.get("bulk_peaks"):
        anno["bulk_peaks"] = pq.read_table(snakemake.input.bulk_peaks).to_pandas()
        anno["bulk_peaks_20kb"] = anno["bulk_peaks"].copy()
        anno["bulk_peaks_20kb"]["Start"] -= 2e4
        anno["bulk_peaks_20kb"]["End"] += 2e4

    rmsk = pr.read_bed(snakemake.input.rmsk, as_df=True)
    rmsk_1kb_3end = pr.read_bed(snakemake.input.rmsk_1kb_3end, as_df=True)
    rmsk_20kb = pr.read_bed(snakemake.input.rmsk_20kb, as_df=True)
    for l1 in ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]:
        anno[l1] = rmsk[rmsk["Name"] == l1]
        anno[f"{l1}_1kb_3end"] = rmsk_1kb_3end[rmsk_1kb_3end["Name"] == l1]
        anno[f"{l1}_20kb"] = rmsk_20kb[rmsk_20kb["Name"] == l1]

    # process windows from each cell of this donor in parallel
    data = Parallel(n_jobs=snakemake.threads)(
        delayed(read_features)(f, anno) for f in snakemake.input.features
    )
    data = pd.concat(data).sort_values(["Chromosome", "Start", "End"])
    data["donor_id"] = snakemake.wildcards.donor

    # save
    data.to_parquet(snakemake.output[0], index=False)

    sys.stderr.close()
