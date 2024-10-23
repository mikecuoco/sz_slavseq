rule make_regions:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        chrom_sizes=config["genome"]["genome"],
    output:
        bed="{outdir}/results/{genome}/{reads}/windows/{donor}/{sample}.bed.gz",
        tabix="{outdir}/results/{genome}/{reads}/windows/{donor}/{sample}.bed.gz.tbi",
    log:
        "{outdir}/results/{genome}/{reads}/windows/{donor}/{sample}.make_regions.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/make_windows.py"


def get_vcfs(wildcards):
    return {
        "xtea": expand(
            rules.vcf2bed.output,
            vcf="xtea",
            allow_missing=True,
        ),
        "megane_percentile": expand(
            rules.vcf2bed.output,
            vcf="megane_percentile",
            allow_missing=True,
        ),
        "megane_gaussian": expand(
            rules.vcf2bed.output,
            vcf="megane_gaussian",
            allow_missing=True,
        ),
        "megane_breakpoints": expand(
            rules.vcf2bed.output,
            vcf="megane_breakpoints",
            allow_missing=True,
        ),
        "graffite": expand(
            rules.vcf2bed.output,
            vcf="graffite",
            allow_missing=True,
        ),
    }


rule bedgraph:
    input:
        bam=rules.sort.output,
        bai=rules.index.output,
    output:
        r1_rr=rules.sort.output[0].replace("bam", "r1_rr.bg"),
        r2_rr=rules.sort.output[0].replace("bam", "r2_rr.bg"),
        contig_rr=rules.sort.output[0].replace("bam", "contig_rr.bg"),
        r1_nonrr=rules.sort.output[0].replace("bam", "r1_nonrr.bg"),
        r2_nonrr=rules.sort.output[0].replace("bam", "r2_nonrr.bg"),
        contig_nonrr=rules.sort.output[0].replace("bam", "contig_nonrr.bg"),
    log:
        rules.sort.log[0].replace("sort", "rr_bedgraph"),
    params:
        min_mapq=20,
    conda:
        "../envs/align.lock.yml"
    shell:
        """
        exec &>> {log}
        samtools view -d RR:1 -F 1024 -f 64 -q {params.min_mapq} -h {input.bam} | genomeCoverageBed -bg -ibam stdin > {output.r1_rr}
        samtools view -d RR:1 -F 1024 -f 128 -q {params.min_mapq} -h {input.bam} | genomeCoverageBed -bg -ibam stdin > {output.r2_rr}
        samtools view -d RR:1 -F 1024 -f 192 -q {params.min_mapq} -h {input.bam} | genomeCoverageBed -bg -ibam stdin > {output.contig_rr}
        samtools view -d RR:0 -F 1024 -f 64 -q {params.min_mapq} -h {input.bam} | genomeCoverageBed -bg -ibam stdin > {output.r1_nonrr}
        samtools view -d RR:0 -F 1024 -f 128 -q {params.min_mapq} -h {input.bam} | genomeCoverageBed -bg -ibam stdin > {output.r2_nonrr}
        samtools view -d RR:0 -F 1024 -f 192 -q {params.min_mapq} -h {input.bam} | genomeCoverageBed -bg -ibam stdin > {output.contig_nonrr}
        """


label_out = {}
for label in [
    "bulk",
    "xtea",
    "megane_gaussian",
    "megane_percentile",
    "megane_breakpoints",
    "graffite",
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
    label_out[label] = rules.make_regions.output.bed.replace(
        ".bed.gz", f".{label}.bed.gz"
    )


# TODO: add closest and overlap labels
rule label:
    input:
        unpack(get_vcfs),
        **rmsk_anno,
        bulk_regions=lambda wc: expand(
            rules.make_regions.output.bed,
            sample=bulk.loc[wc.donor, "sample_id"],
            allow_missing=True,
        ),
        cell_regions=rules.make_regions.output.bed,
        primer_sites=rules.blast_primers.output.bed,
        fai=config["genome"]["fai"],
    output:
        **label_out,
    conda:
        "../envs/ref.lock.yml"
    log:
        rules.make_regions.log[0].replace("make_regions", "label"),
    shell:
        """
        exec &>> {log}

        bedtools intersect -a {input.bulk_regions} -b {input.cell_regions} -wao | bgzip -ci -I {output.bulk}.gzi > {output.bulk}
        bedtools intersect -a {input.primer_sites} -b {input.cell_regions} -wao | bgzip -ci -I {output.primer_sites}.gzi > {output.primer_sites}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.xtea} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.xtea}.gzi > {output.xtea}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_gaussian} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.megane_gaussian}.gzi > {output.megane_gaussian}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_percentile} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.megane_percentile}.gzi > {output.megane_percentile}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_breakpoints} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.megane_breakpoints}.gzi > {output.megane_breakpoints}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.graffite} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.graffite}.gzi > {output.graffite}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1hs} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.l1hs}.gzi > {output.l1hs}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa2} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa2}.gzi > {output.l1pa2}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa3} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa3}.gzi > {output.l1pa3}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa4} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa4}.gzi > {output.l1pa4}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa5} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa5}.gzi > {output.l1pa5}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa6} | bedtools intersect -a - -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa6}.gzi > {output.l1pa6}
        bedtools intersect -a {input.polyA} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa6}.gzi > {output.polyA}
        bedtools intersect -a {input.polyT} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa6}.gzi > {output.polyT}
        """


rule merge_labels:
    input:
        regions=rules.make_regions.output.bed,
        annotations=expand(
            rules.label.output.bulk.replace("bulk", "{annotation}"),
            annotation=label_out,
            allow_missing=True,
        ),
    output:
        rules.label.output.bulk.replace("bulk.bed.gz", "labelled.bed.gz"),
    log:
        rules.label.log[0].replace("label", "merge"),
    conda:
        "../envs/features.yml"
    script:
        "../scripts/merge_labels.py"


rule germline_report:
    input:
        bulk=expand(
            rules.merge_labels.output,
            zip,
            sample=bulk["sample_id"],
            donor=bulk["donor_id"],
            allow_missing=True,
        ),
        megane=expand(
            rules.vcf2bed.output,
            vcf=[
                "xtea",
                "megane_gaussian",
                "megane_percentile",
                "megane_breakpoints",
                "graffite",
            ],
            donor=donors["donor_id"].unique(),
            allow_missing=True,
        ),
        flagstat=expand(
            rules.flagstat.output,
            zip,
            sample=bulk["sample_id"],
            donor=bulk["donor_id"],
            allow_missing=True,
        ),
        meta=config["donors"],
    output:
        rules.reads_report.output[0].replace(
            "reads_report.ipynb", "windows/bulk.bed.gz"
        ),
    conda:
        "../envs/model.yml"
    log:
        notebook=rules.reads_report.output[0].replace(
            "reads_report.ipynb", "windows/germline_report.ipynb"
        ),
    notebook:
        "../scripts/germline_report.py.ipynb"


rule multiinter:
    input:
        bed=expand(
            rules.make_regions.output.bed,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
    output:
        rules.germline_report.output[0].replace("bulk.bed.gz", "multiinter.bed.gz"),
    log:
        rules.germline_report.output[0].replace("bulk.bed.gz", "multiinter.log"),
    conda:
        "../envs/ref.lock.yml"
    shell:
        """
        names=()
        for bed in {input.bed}; do
            names+=($(basename $bed .bed))
        done
        echo ${{names[@]}} > {log}
        bedtools multiinter -i {input.bed} -names ${{names[@]}} > {output} 2>> {log}
        """


rule regions_report:
    input:
        cell_peaks=expand(
            rules.merge_labels.output,
            zip,
            donor=cells["donor_id"],
            sample=cells["sample_id"],
            allow_missing=True,
        ),
        bulk_peaks=rules.germline_report.output,
        multiinter=rules.multiinter.output,
    output:
        rules.germline_report.output[0].replace("bulk.bed.gz", "data.bed.gz"),
    conda:
        "../envs/model.yml"
    log:
        notebook=rules.germline_report.log.notebook.replace(
            "germline_report", "regions_report"
        ),
    notebook:
        "../scripts/regions_report.py.ipynb"
