import sys, pysam

sys.stderr = open(snakemake.log[0], "w")

# iterate over reads
inbam = pysam.AlignmentFile(snakemake.input.bam, "rb")  # must be coordinate sorted
outbam = pysam.AlignmentFile(snakemake.output[0], "wb", template=inbam)
start, end = 0,0
print("Filtering reads...")
for r in inbam.fetch():

	# ignore unmapped reads
	if r.is_unmapped:
		continue

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
		# skip duplicates
		if r.is_reverse:
			if (start == r.reference_end) and (end == r.reference_start) and (r.infer_query_length() == length):
				continue
		else:
			if (start == r.reference_start) and (end == r.reference_end) and (r.infer_query_length() == length):
				continue

		start = r.reference_end if r.is_reverse else r.reference_start
		end = r.reference_start if r.is_reverse else r.reference_end
		length = r.infer_query_length()
		outbam.write(r)

inbam.close()
outbam.close()

sys.stderr.close()
