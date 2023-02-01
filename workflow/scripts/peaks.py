import sys, pysam
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform


sys.stderr = open(snakemake.log[0], "w")

# # iterate over reads
bam = pysam.AlignmentFile(snakemake.input.bam, "rb")  # must be coordinate sorted
peaks = {}
print("Finding L1 reads with same start position...")
for r in bam.fetch():
    # ignore marked duplicates 
    if r.is_duplicate:  
        continue

    # only use primary alignments
    if r.is_secondary or r.is_supplementary:
        continue

    # ignore reads on edges of chrs
    if r.reference_end is None or r.reference_start is None: 
        continue

    if r.reference_name not in peaks.keys():
        peaks[r.reference_name] = {}

    # TODO: check if reads have been clipped
    # TODO: add mate information when adding read2 to peak
    if r.is_read2:  # get read that targets L1

        # get start of read
        if r.is_reverse:
            coord = r.reference_end
        else:
            coord = r.reference_start

        # create/add to peak
        if coord in peaks[r.reference_name].keys():
            peaks[r.reference_name][coord].append(r)
        else:
            peaks[r.reference_name][coord] = [r]


# cluster peaks
print("Clustering pileups...")
t = snakemake.params.cluster_threshold; dist = {}; new_peaks = {}
for chr in peaks.keys():
    
    # skip chr if no peaks
    chr_coords = list(peaks[chr].keys())
    if chr_coords == []:
        continue

    # get distance matrix
    m = np.array([coord for coord in chr_coords]).reshape(-1, 1)
    dist[chr] = squareform(pdist(m, 'euclidean'))

    # cluster peaks
    i = 0; new_peaks[chr] = {} # initialize
    while i < len(chr_coords):
        # get all peaks within t bp
        k,  = np.where(dist[chr][i] <= t)
        coords = [chr_coords[j] for j in k]
        start, end = min(coords), max(coords)
        end += 1 if start == end else 0  # add 1 to end if start == end
        new_peaks[chr][(start, end)] = [r for c in coords for r in peaks[chr][c]]
        i = max(k) + 1


# convert to bed file
print("Saving peaks to BED...")
cols = {"chr": [], "start": [], "end": [], "num_reads": []}
for chr in new_peaks.keys():
    for start, end in new_peaks[chr].keys():
        # skip peaks with only one read
        if len(new_peaks[chr][(start,end)]) == 1:
            continue

        cols["chr"].append(chr)
        cols["start"].append(start)
        cols["end"].append(end)
        cols["num_reads"].append(len(new_peaks[chr][(start,end)]))

# save to file
pd.DataFrame(cols).to_csv(snakemake.output[0], sep="\t", index=False, header=False)

sys.stderr.close()
