#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import sys
from pathlib import Path
import pandas as pd
import pyranges as pr
from pyslavseq.preprocessing import read_cell_features, label
from joblib import Parallel, delayed

sys.stderr = open(snakemake.log[0], "w")

# read blacklist regions
print("Reading blacklist regions", file=sys.stderr)
region_urls = {
    "mhc": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_MHC.bed.gz",
    "kir": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_KIR.bed.gz",
    "trs": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeats_201to10000bp_slop5.bed.gz",
    "segdups": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/SegmentalDuplications/GRCh38_segdups.bed.gz",
    "gaps": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_gaps_slop15kb.bed.gz",
    "false_dup": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_false_duplications_correct_copy.bed.gz",
}

regions = {}
for id, url in region_urls.items():
    regions[id] = pd.read_csv(
        url,
        sep="\t",
        header=None,
        skiprows=1,
        names=["Chromosome", "Start", "End"],
    )
blacklist = pd.concat([r for _, r in regions.items()])

# read annotations
print("Reading L1 annotations", file=sys.stderr)
anno = {
    "xtea": pr.read_bed(snakemake.input.xtea).df,
    "xtea_1kb_3end": pr.read_bed(snakemake.input.xtea_1kb_3end).df,
    "xtea_20kb": pr.read_bed(snakemake.input.xtea_20kb).df,
}

rmsk = pr.read_bed(snakemake.input.rmsk).df
rmsk_1kb_3end = pr.read_bed(snakemake.input.rmsk_1kb_3end).df
rmsk_20kb = pr.read_bed(snakemake.input.rmsk_20kb).df

for l1 in ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]:
    anno[l1] = rmsk.loc[rmsk["Name"].str.contains(l1), :]
    anno[l1 + "_1kb_3end"] = rmsk_1kb_3end.loc[
        rmsk_1kb_3end["Name"].str.contains(l1), :
    ]
    anno[l1 + "_20kb"] = rmsk_20kb.loc[rmsk_20kb["Name"].str.contains(l1), :]


def read_data(files, n_jobs: int):
    """Reads data from each cell in parallel"""
    # check log files of features to make sure they finished
    # for f in files:
    #     log = Path(f.rstrip("_windows.pqt")).with_suffix(".log")
    #     with open(log, "r") as f:
    #         # assert last line contains "Done"
    #         assert "Done" in f.readlines()[-1]

    data = Parallel(n_jobs=n_jobs)(
        delayed(read_cell_features)(
            f, cell_id=Path(f).stem.rstrip("_windows"), exclude_ref=False
        )
        for f in files
    )
    data = pd.concat(data).sort_values(["Chromosome", "Start", "End"])

    print("Removing non-autosomal regions", file=sys.stderr)
    data = data.loc[data["Chromosome"].isin([f"chr{i}" for i in range(1, 23)])]

    # remove blacklist regions
    print("Removing blacklist regions", file=sys.stderr)
    data = label(data, blacklist, "blacklist")
    data = data.loc[data["blacklist"] == False]
    data = data.drop(columns=["blacklist"])  # remove blacklist column

    return data


# LABEL WINDOWS
print("Reading windows", file=sys.stderr)
data = read_data(snakemake.input.windows, snakemake.threads)

# label knrgls
print("Labeling L1 regions in windows", file=sys.stderr)
# NOT WORKING: add bulk peaks if available
# if snakemake.input.get("bulk_peaks"):
# 	anno["bulk_peaks"] = pq.read_table(snakemake.input.bulk_peaks).to_pandas()
# 	anno["bulk_peaks_20kb"] = anno["bulk_peaks"].copy()
# 	anno["bulk_peaks_20kb"]["End"] += 2e4

for id, df in anno.items():
    data = label(data, df, id)
    assert (
        data.shape[0]
        == data[["Chromosome", "Start", "End", "cell_id"]].drop_duplicates().shape[0]
    ), f"some rows have been duplicated during {id} labeling!"

data["donor_id"] = snakemake.wildcards.donor

# add motif scores
en_motif = pd.read_parquet(snakemake.input.en_motif).set_index(
    ["Chromosome", "Start", "End"]
)
data = (
    data.set_index(["Chromosome", "Start", "End"])
    .join(en_motif, how="left")
    .reset_index()
)

# save
print("Saving labelled windows", file=sys.stderr)
data.to_parquet(snakemake.output.windows, index=False)
data[data["n_ref_reads"] == 0].to_parquet(
    snakemake.output.windows_nonrefonly, index=False
)
del data

# LABEL PEAKS
print("Reading peaks", file=sys.stderr)
data = read_data(snakemake.input.peaks, snakemake.threads)

# label knrgls
print("Labeling L1 regions in peaks", file=sys.stderr)
for id, df in anno.items():
    if id == "xtea":
        data = label(data, df, id, add_id=True, slack=50)
    else:
        data = label(data, df, id)
    assert (
        data.shape[0]
        == data[["Chromosome", "Start", "End", "cell_id"]].drop_duplicates().shape[0]
    ), f"some rows have been duplicated during {id} labeling!"

data["donor_id"] = snakemake.wildcards.donor

# save
print("Saving labelled peaks", file=sys.stderr)
data.to_parquet(snakemake.output.peaks, index=False)
data[data["n_ref_reads"] == 0].to_parquet(
    snakemake.output.peaks_nonrefonly, index=False
)

sys.stderr.close()
