import sys, pysam
import pandas as pd
import numpy as np

sys.stderr = open(snakemake.log[0], "w")

# iterate over reads
bam = pysam.AlignmentFile(snakemake.input.bam, "rb")  # must be coordinate sorted
peaks = {}
print("Calling non-reference L1 insertions...")
for r in bam.fetch():

    # ignore unmapped reads
    if r.is_unmapped:
        continue

    # make peaks from read1 if
    # 1. is not a proper pair
    # 2. is discordant (template length <= 0, or > 2000)
    # 3. has clipped bases (S or H in MC (mate-cigar) tag)
    if r.is_read1:
        if not r.is_proper_pair or r.template_length <= 0 or r.template_length > 2000 or "S" in r.get_tag("MC") or "H" in r.get_tag("MC"):

            # if overlap with previous peak, add to peak
            if r.reference_name in peaks.keys() and r.get_overlap(start, end) > 0: 
                peaks[r.reference_name][-1].append(r)
                if r.is_reverse:
                    start = r.reference_end if r.reference_end < start else start
                    end = r.reference_start if r.reference_start > end else end
                else: 
                    start = r.reference_start if r.reference_start < start else start
                    end = r.reference_end if r.reference_end > end else end
                continue

            # create new list for each chr
            if r.reference_name not in peaks.keys():
                peaks[r.reference_name] = []
                
            # if no overlap, create new peak
            peaks[r.reference_name].append([r])
            start = r.reference_end if r.is_reverse else r.reference_start
            end = r.reference_start if r.is_reverse else r.reference_end

# convert to bed file
print("Saving peaks to BED...")
bed = {"chr": [], "start": [], "end": [], "r1_count": [], "r1_uniq_starts": [], "r1_uniq_ends": [], "r1_primary": [], "r1_supplementary": []}
for chr in peaks.keys():
    for p in peaks[chr]:
        # skip peaks with only one read
        if len(p) == 1:
            continue

        features = {"r1_start": np.array([], dtype=int), "r1_end": np.array([], dtype=int), "r1_primary": 0, "r1_supplementary": 0}
        for r in p:
            features["r1_start"] = np.append(features["r1_start"], r.reference_end if r.is_reverse else r.reference_start)
            features["r1_end"] = np.append(features["r1_end"], r.reference_start if r.is_reverse else r.reference_end)
            
            if r.is_supplementary:
                features["r1_supplementary"] += 1
            else:
                features["r1_primary"] += 1

        bed["chr"].append(chr)
        bed["start"].append(features["r1_start"].min())
        bed["end"].append(features["r1_end"].max())
        bed["r1_count"].append(len(p))
        bed["r1_uniq_starts"].append(len(np.unique(features["r1_start"])))
        bed["r1_uniq_ends"].append(len(np.unique(features["r1_end"])))
        bed["r1_primary"].append(features["r1_primary"])
        bed["r1_supplementary"].append(features["r1_supplementary"])

# save to file
pd.DataFrame(bed).to_csv(snakemake.output[0], sep="\t", index=False, header=False)

sys.stderr.close()
