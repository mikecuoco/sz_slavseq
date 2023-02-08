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
    # 1. is read1
    # 2. has YA (R2 L1 alignment score) and YG (R2 reference genome alignment score) tags
    # 3. YA > YG
    if (
        r.is_read1
        and r.has_tag("YA")
        and r.has_tag("YG")
        and r.get_tag("YA") > r.get_tag("YG") # ignore reference L1 insertions
        and r.get_tag("YA") > 20
        and not r.has_tag("SA")
        and not r.has_tag("XA")
        and r.mapping_quality >= 60
    ):
        

        # if overlap with previous peak, add to peak
        if (r.reference_name in peaks.keys()) and (r.get_overlap(start, end) > 10):

            # skip duplicates
            if r.is_reverse and (r.reference_end == start) and (r.reference_start == end):
                continue
            elif (not r.is_reverse) and (r.reference_start == start) and (r.reference_end == end):
                continue


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

        # create new peak
        peaks[r.reference_name].append([r])
        start = r.reference_end if r.is_reverse else r.reference_start
        end = r.reference_start if r.is_reverse else r.reference_end
bam.close()

# convert to bed file
print("Collecting peak features...")
# TODO: iterate over peaks and collect features from all reads in each peak
bed = {
    "chr": [],
    "start": [],
    "end": [],
    "width": [],
    "total_reads": [],
}

# reopen bam file
bam = pysam.AlignmentFile(snakemake.input.bam, "rb")  # must be coordinate sorted
for chr in peaks.keys():
    for p in peaks[chr]:
        # skip peaks with only one read
        if len(p) == 1:
            continue

        features = {
            "r1_start": np.array([], dtype=int),
            "r1_end": np.array([], dtype=int),
        }
        for r in p:
            features["r1_start"] = np.append(
                features["r1_start"],
                r.reference_end if r.is_reverse else r.reference_start,
            )
            features["r1_end"] = np.append(
                features["r1_end"],
                r.reference_start if r.is_reverse else r.reference_end,
            )

        bed["chr"].append(chr)
        bed["start"].append(features["r1_start"].min())
        bed["end"].append(features["r1_end"].max())
        bed["width"].append(bed["end"][-1] - bed["start"][-1])
        bed["total_reads"].append(len(p))

        # # find peak summit
        # covr = np.zeros(bed["width"][-1])
        # for c in bam.count_coverage(chr, bed["start"][-1], bed["end"][-1]):
        #     covr = np.add(covr,c)
        # summit = np.where(covr == max(covr))
        # import pdb; pdb.set_trace()
        # bed["start"][-1] = bed["start"][-1] + summit[0][0]
        # bed["end"][-1] = bed["start"][-1] + summit[0][-1]
        # bed["width"][-1] = bed["end"][-1] - bed["start"][-1]

# save to file
print("saving to file...")
pd.DataFrame(bed).to_csv(snakemake.output[0], index=False, sep="\t", header=None)

sys.stderr.close()
