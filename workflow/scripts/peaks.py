import sys, pysam
import pandas as pd


sys.stderr = open(snakemake.log[0], "w")

# iterate over reads
bam = pysam.AlignmentFile(snakemake.input.bam, "rb")  # must be coordinate sorted
peaks = {}
print("Calling non-reference L1 insertions...")
for r in bam.fetch():

    # # only use primary alignments
    # if r.is_secondary or r.is_supplementary or r.is_unmapped:
    #     continue

    # ignore unmapped reads
    if r.is_unmapped:
        continue

    # investigate discordant read1
    if r.is_read1:
        if r.template_length == 0 or abs(r.template_length) > 2000:
            if r.reference_name not in peaks.keys():  # create new list for each chr
                peaks[r.reference_name] = []
                i, start, end = 0, 0, 0
                peaks[r.reference_name].append([r])
            elif (
                r.get_overlap(start, end) > 0
            ):  # if overlaps last read, add to existing peak
                peaks[r.reference_name][i].append(r)
            else:  # create new peak
                i += 1
                peaks[r.reference_name].append([r])
            start = r.reference_end if r.is_reverse else r.reference_start
            end = r.reference_start if r.is_reverse else r.reference_end


# convert to bed file
print("Saving peaks to BED...")
bed = {"chr": [], "start": [], "end": [], "num_reads": []}
for chr in peaks.keys():
    for p in peaks[chr]:
        # skip peaks with only one read
        if len(p) == 1:
            continue

        bed["chr"].append(chr)
        bed["start"].append(
            min([r.reference_end if r.is_reverse else r.reference_start for r in p])
        )
        bed["end"].append(
            max([r.reference_start if r.is_reverse else r.reference_end for r in p])
        )
        bed["num_reads"].append(len(p))

# save to file
pd.DataFrame(bed).to_csv(snakemake.output[0], sep="\t", index=False, header=False)

sys.stderr.close()
