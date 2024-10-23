# rule bigwig:
# 	input:
# 		rules.sort.output,
# 	output:
# 		r1=rules.sort.output[0].replace("bam", "r1.bw"),
# 		r2=rules.sort.output[0].replace("bam", "r2.bw"),
# 		contig=rules.sort.output[0].replace("bam", "contig.bw"),
# 	log:
# 		rules.sort.log[0].replace("sort", "bigwig"),
# 	params:
# 		extra=f"""
# 		--minMappingQuality 20 \
# 		--effectiveGenomeSize {config['genome']['size']} \
# 		--ignoreDuplicates \
# 		--normalizeUsing RPGC \
# 		--exactScaling \
# 		""",
# 	conda:
# 		"../envs/align.lock.yml"
# 	shell:
# 		"""
# 		bamCoverage -b {input} -o {output.r1} -p {threads} {params.extra} --samFlagInclude 64 2> {log}
# 		bamCoverage -b {input} -o {output.r2} -p {threads} {params.extra} --samFlagInclude 128 2> {log}
# 		bamCoverage -b {input} -o {output.contig} -p {threads} {params.extra} --samFlagInclude 192 2> {log}
# 		"""
