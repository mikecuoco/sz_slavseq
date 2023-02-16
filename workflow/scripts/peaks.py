import sys, pysam
import pandas as pd
import numpy as np

sys.stderr = open(snakemake.log[0], "w")


class Peak:
    """Store collection of reads in a peak"""

    def __init__(self, r):
        self.chr = r.reference_name
        self.start = r.reference_start
        self.end = r.reference_end
        self.is_reverse = not r.is_reverse if r.is_read1 else r.is_reverse
        self.width = self.end - self.start
        if r.is_read1:
            self.r1 = [r]
            self.r1_ends = [r.query_alignment_end]
            self.r2 = []
            self.r2_starts = []
        elif r.is_read2:
            self.r1 = []
            self.r1_ends = []
            self.r2 = [r]
            self.r2_starts = [r.query_alignment_start]
        self.coverage = np.ones(r.reference_end - r.reference_start)
        assert self.end >= self.start

    def add_read(self, r):
        """Add read to peak, extending peak if necessary"""

        if r.is_read1:
            self.r1.append(r)
            self.r1_ends.append(r.query_alignment_end)
        elif r.is_read2:
            self.r2.append(r)
            self.r2_starts.append(r.query_alignment_start)

        # extend peak
        if r.reference_end > self.end:
            self.coverage = np.append(
                self.coverage, np.zeros(r.reference_end - self.end)
            )
            self.end = r.reference_end
            self.width = self.end - self.start

        start = r.reference_start - self.start
        end = r.reference_end - self.start
        self.coverage[start:end] += 1

        assert self.end >= self.start

    # NOTE: this is not used in the current version of the pipeline,
    #       it has the potential to remove reads that are not true duplicates
    def is_duplicate(self, r):
        """Check if read is duplicate of previous read"""
        if (
            r.is_read1
            and len(self.r1) > 0
            and (r.reference_name == self.chr)
            and (r.query_alignment_start == self.r1[-1].query_alignment_start)
            and (r.query_alignment_end == self.r1[-1].query_alignment_end)
        ):
            return True
        elif (
            r.is_read2
            and len(self.r2) > 0
            and (r.reference_name == self.chr)
            and (r.query_alignment_start == self.r2[-1].query_alignment_start)
            and (r.query_alignment_end == self.r2[-1].query_alignment_end)
        ):
            return True
        else:
            return False

    def has_mate(self, r):
        """Check if read has a mate in the peak"""
        if r.is_read2:
            for r1 in self.r1:
                if r.query_name == r1.query_name:
                    return True
        elif r.is_read1:
            for r2 in self.r2:
                if r.query_name == r2.query_name:
                    return True


# iterate over reads
bam = pysam.AlignmentFile(snakemake.input.bam, "rb")  # must be coordinate sorted
bed = {
    "chr": [],
    "start": [],
    "end": [],
    "width": [],
    "nreads": [],
}
p = None
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
        and (r.has_tag("YA") and r.has_tag("YG"))
        and (r.get_tag("YA") > 20 and r.get_tag("YA") > r.get_tag("YG"))
        and (not r.has_tag("XA"))
        and (r.mapping_quality >= 60)
    ):

        # initialize
        if p is None:
            p = Peak(r)

        # if has mate in peak or overlap with peak, add to peak
        elif (
            (r.reference_name == p.chr)
            and (
                (r.is_read1 and (r.is_reverse != p.is_reverse))
                or (r.is_read2 and (r.is_reverse == p.is_reverse))
            )
            and (r.get_overlap(p.start, p.end) > 0 or p.has_mate(r))
        ):
            p.add_read(r)

        # else, make new peak
        else:
            # add peak to bed if more than 1 read
            if len(p.r1) > 1:
                bed["chr"].append(p.chr)
                bed["start"].append(p.start)
                bed["end"].append(p.end)
                bed["width"].append(p.width)
                bed["nreads"].append(len(p.r1) + len(p.r2))
            p = Peak(r)

bam.close()

# save to file
print("saving to file...")
pd.DataFrame(bed).to_csv(snakemake.output.peaks, index=False, sep="\t", header=None)

sys.stderr.close()
