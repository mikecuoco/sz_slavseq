#!/usr/bin/env python
__author__ = "Rohini Gadde", "Michael Cuoco"

import sys
from pathlib import Path
import pyarrow.parquet as pq
import pandas as pd
import pyranges as pr
from pyslavseq.preprocessing import read_cell_features, label
from joblib import Parallel, delayed


sys.stderr = open(snakemake.log[0], "w")

# check log files of features to make sure they finished
for f in snakemake.input.features:
    log = Path(f.rstrip("_windows.pqt")).with_suffix(".log")
    with open(log, "r") as f:
        # assert last line contains "Done"
        assert "Done" in f.readlines()[-1]

# process windows from each cell of this donor in parallel
data = Parallel(n_jobs=snakemake.threads)(
    delayed(read_cell_features)(
        f, cell_id=Path(f).stem.rstrip("_windows"), exclude_ref=False
    )
    for f in snakemake.input.features
)
data = pd.concat(data).sort_values(["Chromosome", "Start", "End"])

# keep autosomes
print("Removing non-autosomal windows", file=sys.stderr)
data = data.loc[data["Chromosome"].isin([f"chr{i}" for i in range(1, 23)])]

# remove blacklist regions
print("Removing blacklist regions", file=sys.stderr)
mhc = pd.read_csv(
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_MHC.bed.gz",
    sep="\t",
    header=None,
    skiprows=1,
    names=["Chromosome", "Start", "End"],
)
kir = pd.read_csv(
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_KIR.bed.gz",
    sep="\t",
    header=None,
    skiprows=1,
    names=["Chromosome", "Start", "End"],
)
trs = pd.read_csv(
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/LowComplexity/GRCh38_AllTandemRepeats_201to10000bp_slop5.bed.gz",
    sep="\t",
    header=None,
    skiprows=1,
    names=["Chromosome", "Start", "End"],
)
segdups = pd.read_csv(
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/SegmentalDuplications/GRCh38_segdups.bed.gz",
    sep="\t",
    header=None,
    skiprows=1,
    names=["Chromosome", "Start", "End"],
)
gaps = pd.read_csv(
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_gaps_slop15kb.bed.gz",
    sep="\t",
    header=None,
    skiprows=1,
    names=["Chromosome", "Start", "End"],
)
false_dup = pd.read_csv(
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.0/GRCh38/OtherDifficult/GRCh38_false_duplications_correct_copy.bed.gz",
    sep="\t",
    header=None,
    skiprows=1,
    names=["Chromosome", "Start", "End"],
)
blacklist = pd.concat([mhc, trs, segdups, gaps, false_dup, kir])
data = label(data, blacklist, "blacklist")
data = data.loc[data["blacklist"] == False]

# remove blacklist column
data = data.drop(columns=["blacklist"])

# label knrgls
print("Labeling knrgls", file=sys.stderr)

anno = {
    "xtea": pr.read_bed(snakemake.input.xtea).df,
    "xtea_1kb_3end": pr.read_bed(snakemake.input.xtea_1kb_3end).df,
    "xtea_20kb": pr.read_bed(snakemake.input.xtea_20kb).df,
}

# add bulk peaks if available
# if snakemake.input.get("bulk_peaks"):
# 	anno["bulk_peaks"] = pq.read_table(snakemake.input.bulk_peaks).to_pandas()
# 	anno["bulk_peaks_20kb"] = anno["bulk_peaks"].copy()
# 	anno["bulk_peaks_20kb"]["End"] += 2e4

rmsk = pr.read_bed(snakemake.input.rmsk).df
rmsk_1kb_3end = pr.read_bed(snakemake.input.rmsk_1kb_3end).df
rmsk_20kb = pr.read_bed(snakemake.input.rmsk_20kb).df

for l1 in ["L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6"]:
    anno[l1] = rmsk.loc[rmsk["Name"].str.contains(l1), :]
    anno[l1 + "_1kb_3end"] = rmsk_1kb_3end.loc[
        rmsk_1kb_3end["Name"].str.contains(l1), :
    ]
    anno[l1 + "_20kb"] = rmsk_20kb.loc[rmsk_20kb["Name"].str.contains(l1), :]

for id, df in anno.items():
    data = label(data, df, id)

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

# check that no rows have been duplicated
assert (
    data.shape[0]
    == data[["Chromosome", "Start", "End", "cell_id"]].drop_duplicates().shape[0]
), "some rows have been duplicated during labeling!"

# save
data.to_parquet(snakemake.output[0], index=False)

# remove ref
data[data["n_ref_reads"] == 0].to_parquet(snakemake.output[1], index=False)

sys.stderr.close()
