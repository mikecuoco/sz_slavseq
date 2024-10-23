rule make_bigwig:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        chrom_sizes=config["genome"]["genome"],
    output:
        bg="{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.bg",
        bw="{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.bw",
    log:
        "{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.make_bigwig.log",
    conda:
        "../envs/ref.lock.yml"
    shell:
        """
        samtools view -h -F 1412 {input.bam} | bedtools genomecov -bg -ibam stdin -3 | sort -k1,1 -k2,2n > {output.bg} 2> {log}
        bedGraphToBigWig {output.bg} {input.chrom_sizes} {output.bw} 2>> {log}
        """


rule merge_donor_bigwig:
    input:
        bw=expand(
            rules.make_bigwig.output.bw,
            sample=lambda wc: cells[cells.donor_id == wc.donor_id]["sample_id"],
            allow_missing=True,
        ),
        chrom_sizes=config["genome"]["genome"],
    output:
        bg=temp("{outdir}/results/{genome}/{reads}/bigwig/{donor}/merged.bg"),
        bg_sorted=temp(
            "{outdir}/results/{genome}/{reads}/bigwig/{donor}/merged.sorted.bg"
        ),
        bw="{outdir}/results/{genome}/{reads}/bigwig/{donor}/merged.bw",
    log:
        "{outdir}/results/{genome}/{reads}/bigwig/{donor}/merge.log",
    conda:
        "../envs/ref.lock.yml"
    shell:
        """
        bigWigMerge -max -threshold=10 {input.bw} {output.bg} 2> {log}
        sort -k1,1 -k2,2bn {output.bg} > {output.bg_sorted} 2>> {log}
        bedGraphToBigWig {output.bg_sorted} {input.chrom_sizes} {output.bw} 2>> {log}
        """


rule merge_bigwig:
    input:
        bw=expand(
            rules.make_bigwig.output.bw,
            zip,
            sample=cells["sample_id"],
            donor=cells["donor_id"],
            allow_missing=True,
        ),
        chrom_sizes=config["genome"]["genome"],
    output:
        bg=temp("{outdir}/results/{genome}/{reads}/bigwig/all.bg"),
        bg_sorted=temp("{outdir}/results/{genome}/{reads}/bigwig/all.sorted.bg"),
        bw="{outdir}/results/{genome}/{reads}/bigwig/all.bw",
    log:
        "{outdir}/results/{genome}/{reads}/bigwig/merge_bigwig.log",
    conda:
        "../envs/ref.lock.yml"
    shell:
        """
        bigWigMerge -max -threshold=10 {input.bw} {output.bg} 2> {log}
        sort -k1,1 -k2,2bn {output.bg} > {output.bg_sorted} 2>> {log}
        bedGraphToBigWig {output.bg_sorted} {input.chrom_sizes} {output.bw} 2>> {log}
        """


rule greedy:
    input:
        bw=rules.merge_bigwig.output.bw,
    output:
        bed="{outdir}/results/{genome}/{reads}/bigwig/peaks.bed.gz",
    log:
        "{outdir}/results/{genome}/{reads}/bigwig/greedy.log",
    conda:
        "../envs/features.yml"
    params:
        blacklist_bp=2500,
        extend_bp=100,
    script:
        "../scripts/greedy_peaks.py"


rule get_features:
    input:
        peaks=rules.greedy.output.bed,
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
    output:
        bed="{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.bed.gz",
        tbi="{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.bed.gz.tbi",
    log:
        "{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.get_features.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_features.py"


rule merge_greedy_windows:
    input:
        expand(
            rules.get_features.output.bed,
            zip,
            sample=cells["sample_id"],
            donor=cells["donor_id"],
            allow_missing=True,
        ),
    output:
        bed="{outdir}/results/{genome}/{reads}/bigwig/data.bed.gz",
        tbi="{outdir}/results/{genome}/{reads}/bigwig/data.bed.gz.tbi",
    conda:
        "../envs/ref.lock.yml"
    log:
        "{outdir}/results/{genome}/{reads}/bigwig/merge.log",
    shell:
        """
        # remove .gz extension from output files
        bed=$(echo {output.bed} | sed 's/.gz//g')

        # run the code
        zcat {input[0]} | grep "^#" > $bed
        xargs -n 1 zcat <<< "{input}" | grep -v "^#" | sort -k1,1 -k2,2n >> $bed
        bgzip -f $bed
        tabix -p bed {output.bed}
        """


label_out = {}
for label in [
    "megane",
    "primer_sites",
    "l1hs",
    "l1pa2",
    "l1pa3",
    "l1pa4",
    "l1pa5",
    "l1pa6",
    "polyA",
    "polyT",
]:
    label_out[label] = rules.greedy.output.bed.replace(".bed.gz", f".{label}.bed.gz")


# TODO: add closest and overlap labels
rule label_greedy_windows:
    input:
        **rmsk_anno,
        **knrgl_anno,
        regions=rules.greedy.output.bed,
        primer_sites=rules.blast_primers.output.bed,
        fai=config["genome"]["fai"],
    output:
        **label_out,
    conda:
        "../envs/ref.lock.yml"
    log:
        rules.greedy.log[0].replace("greedy", "label"),
    shell:
        """
        exec &>> {log}

        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.primer_sites} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.primer_sites}.tbi > {output.primer_sites}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.megane}.tbi > {output.megane}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1hs} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.l1hs}.tbi > {output.l1hs}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa2} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.l1pa2}.tbi > {output.l1pa2}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa3} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.l1pa3}.tbi > {output.l1pa3}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa4} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.l1pa4}.tbi > {output.l1pa4}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa5} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.l1pa5}.tbi > {output.l1pa5}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa6} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.l1pa6}.tbi > {output.l1pa6}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.polyA} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.polyA}.tbi > {output.polyA}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.polyT} | bedtools intersect -a - -b {input.regions} -wao | bgzip -ci -I {output.polyT}.tbi > {output.polyT}
        """


rule merge_greedy_window_labels:
    input:
        regions=rules.greedy.output.bed,
        annotations=expand(
            rules.label_greedy_windows.output.megane.replace("megane", "{annotation}"),
            annotation=label_out,
            allow_missing=True,
        ),
    output:
        bed=rules.label_greedy_windows.output.megane.replace(
            "megane.bed.gz", "labelled.bed.gz"
        ),
    log:
        rules.label_greedy_windows.log[0].replace("label", "merge_labels"),
    conda:
        "../envs/features.yml"
    script:
        "../scripts/merge_labels.py"


# rule greedy_windows_report:
# 	input:
# 		labels=rules.merge_greedy_window_labels.output.bed,
# 		cells=rules.merge_greedy_windows.output.bed,
# 		bulk=rules.bulk_peaks_report.output.bed,
# 		meta=config["donors"],
# 	output:
# 		"{outdir}/results/{genome}/{reads}/bigwig/data_final.bed.gz",
# 	conda:
# 		"../envs/model.yml"
# 	log:
# 		notebook="{outdir}/results/{genome}/{reads}/bigwig/report.ipynb",
# 	notebook:
# 		"../scripts/greedy_windows_report.py.ipynb"


# rule model_greedy_windows:
# 	input:
# 		data=rules.greedy_windows_report.output,
# 	output:
# 		model=rules.greedy_windows_report.output[0].replace("data.bed.gz", "model.pkl"),
# 		best_hp=rules.greedy_windows_report.output[0].replace(
# 			"data.bed.gz", "best_hp.json"
# 		),
# 		history=rules.greedy_windows_report.output[0].replace(
# 			"data.bed.gz", "history.log"
# 		),
# 		predictions=rules.greedy_windows_report.output[0].replace(
# 			"data.bed.gz", "predictions.pqt"
# 		),
# 		shap_values=rules.greedy_windows_report.output[0].replace(
# 			"data.bed.gz", "shap_values.pqt"
# 		),
# 	conda:
# 		"../envs/model_gpu.yml" if shutil.which("nvidia-smi") else "../envs/model.yml"
# 	params:
# 		max_iter=100,
# 		random_state=1,
# 		# features=lambda wc: feature_sets[wc.features],
# 	threads: 1
# 	log:
# 		notebook=rules.greedy_windows_report.log.notebook.replace(
# 			"report", "tune"
# 		),
# 	notebook:
# 		"../scripts/tune.py.ipynb"


rule greedy_windows:
    input:
        expand(
            rules.merge_greedy_windows.output,
            reads="filtered",
            genome=config["genome"]["name"],
            outdir=config["outdir"],
        ),
