import sys, pysam
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


class Peak:
    """Store collection of reads in a peak"""

    def __init__(self, r):
        self.chr = r.reference_name
        self.start = r.reference_start
        self.end = r.reference_end
        self.is_reverse = r.is_reverse
        self.reads = [r]
        assert self.end >= self.start

    def add_read(self, r):
        self.reads.append(r)
        self.end = r.reference_end if r.reference_end > self.end else self.end
        assert self.end >= self.start

    def is_duplicate(self, r):
        if (
            (r.reference_name == self.chr)
            and (r.reference_start == self.reads[-1].reference_start)
            and (r.reference_end == self.reads[-1].reference_end)
        ):
            return True
        else:
            return False


# iterate over reads
bam = pysam.AlignmentFile(snakemake.input.bam, "rb")  # must be coordinate sorted
peaks = []
print("Calling non-reference L1 insertion peaks...")
for r in bam.fetch():

    # make peaks from read1 if
    # 1. is a read1 primary alignment
    # 2. mate aligns better to L1 than reference genome (YA >20, YG < 10)
    # 3. does not have any secondary alignments (SA, XA)
    # 4. has high mapping quality (60)
    if (
        r.is_read1
        and (not r.is_unmapped)
        and (not r.is_secondary)
        and (not r.is_supplementary)
        and r.has_tag("YA")
        and r.has_tag("YG")
        and (r.get_tag("YA") > 20)
        and (r.get_tag("YG") < 10)
        and (not r.has_tag("SA"))
        and (not r.has_tag("XA"))
        and (r.mapping_quality >= 60)
    ):
        # initialize
        if len(peaks) == 0:
            peaks.append(Peak(r))

        # if overlap with previous peak, add to peak
        elif (
            (r.reference_name == peaks[-1].chr)
            and (r.is_reverse == peaks[-1].is_reverse)
            and (r.get_overlap(peaks[-1].start, peaks[-1].end) > 0)
            # and (abs(r.reference_end - peaks[-1].reads[-1].reference_start) < 50)
        ):
            # skip duplicates
            if not peaks[-1].is_duplicate(r):
                peaks[-1].add_read(r)

        # else, make new peak
        else:
            # import pdb; pdb.set_trace()
            peaks.append(Peak(r))

bam.close()

# convert to bed file
print("Defining peak boundaries...")

bed = {
    "chr": [],
    "start": [],
    "end": [],
    "width": [],
    "total_reads": [],
}

for p in peaks:
    # skip peaks with only one read
    if len(p.reads) == 1:
        continue

    bed["chr"].append(p.chr)
    bed["start"].append(p.start)
    bed["end"].append(p.end)
    bed["width"].append(bed["end"][-1] - bed["start"][-1])
    bed["total_reads"].append(len(p.reads))

# save to file
print("saving to file...")
pd.DataFrame(bed).to_csv(snakemake.output[0], index=False, sep="\t", header=None)

sys.stderr.close()
