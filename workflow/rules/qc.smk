def get_fastqc_input(wc):
    if wc.fastq == "raw":
        reads = [
            samples.loc[wc.donor, wc.sample][["R1"]],
            samples.loc[wc.donor, wc.sample][["R2"]],
        ]
    elif wc.fastq == "trimmed":
        reads = [rules.trim_adapters.output.fastq1, rules.trim_adapters.output.fastq2]
    elif wc.fastq == "filtered":
        reads = [rules.filter_read2.output.fastq1, rules.filter_read2.output.fastq2]
    return reads[0] if wc.read == "R1" else reads[1]


rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="{outdir}/results/qc/fastqc/{donor}/{sample}_{read}.{fastq}.html",
        zip="{outdir}/results/qc/fastqc/{donor}/{sample}_{read}.{fastq}.fastqc.zip",
    log:
        "{outdir}/results/qc/fastqc/{donor}/{sample}_{read}.{fastq}.log",
    wrapper:
        "v1.21.0/bio/fastqc"


rule flagstat:
    input:
        "{outdir}/results/{genome}/align/{donor}/{sample}.{stage}.bam",
    output:
        "{outdir}/results/{genome}/qc/flagstat/{donor}/{sample}.{stage}.flagstat",
    log:
        "{outdir}/results/{genome}/qc/flagstat/{donor}/{sample}.{stage}.flagstat.log",
    wrapper:
        "v1.21.0/bio/samtools/flagstat"


def get_l1_coverage_anno(wildcards):
    if "rmsk" in wildcards.anno:
        return rules.rmsk_to_bed.output[wildcards.anno]
    elif "xtea" in wildcards.anno:
        return rules.xtea_to_bed.output[wildcards.anno]
    elif "bulk_peaks" in wildcards.anno:
        return rules.call_bulk_peaks.output[wildcards.anno]
    else:
        raise ValueError("Unknown annotation: {}".format(wildcards.anno))


rule l1_coverage:
    input:
        bam=rules.sambamba_sort.output,
        bai=rules.sambamba_index.output,
        anno=get_l1_coverage_anno,
    output:
        r1="{outdir}/results/{genome}/qc/l1_coverage/{donor}/{sample}.{anno}.r1.txt",
        r2="{outdir}/results/{genome}/qc/l1_coverage/{donor}/{sample}.{anno}.r2.txt",
    log:
        "{outdir}/results/{genome}/qc/l1_coverage/{donor}/{sample}.{anno}.log",
    params:
        extra=" ",  # optional additional parameters as string
    conda:
        "../envs/ref.yml"
    shell:
        """
        samtools view -b -f 64 {input.bam} | bedtools coverage -a {input.anno} -b stdin > {output.r1} 2> {log}
        samtools view -b -f 128 {input.bam} | bedtools coverage -a {input.anno} -b stdin > {output.r2} 2>> {log}
        """


rule qc_summary:
    input:
        rmsk_1kb_3end=rules.rmsk_to_bed.output.rmsk_1kb_3end,
        trim_qc=expand(
            rules.trim_adapters.output.qc,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
        filter_qc=expand(
            rules.filter_read2.output.qc,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
        fastqc_zip=expand(
            expand(
                rules.fastqc.output.zip,
                zip,
                sample=samples["sample_id"],
                donor=samples["donor_id"],
                allow_missing=True,
            ),
            fastq=["trimmed"],
            read=["R1", "R2"],
            allow_missing=True,
        ),
        flagstat=expand(
            expand(
                rules.flagstat.output,
                zip,
                sample=samples["sample_id"],
                donor=samples["donor_id"],
                allow_missing=True,
            ),
            ref=["genome"],
            allow_missing=True,
        ),
        l1_coverage_r1=expand(
            expand(
                rules.l1_coverage.output.r1,
                zip,
                sample=samples["sample_id"],
                donor=samples["donor_id"],
                allow_missing=True,
            ),
            anno=[
                "knrgl",
                "knrgl_1kb_3end",
                "knrgl_20kb",
                "rmsk",
                "rmsk_1kb_3end",
                "rmsk_20kb",
            ],
            allow_missing=True,
        ),
        l1_coverage_r2=expand(
            expand(
                rules.l1_coverage.output.r2,
                zip,
                sample=samples["sample_id"],
                donor=samples["donor_id"],
                allow_missing=True,
            ),
            anno=[
                "knrgl",
                "knrgl_1kb_3end",
                "knrgl_20kb",
                "rmsk",
                "rmsk_1kb_3end",
                "rmsk_20kb",
            ],
            allow_missing=True,
        ),
        bams=expand(
            rules.sambamba_sort.output,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/qc/qc_summary.ipynb",
    log:
        "{outdir}/results/qc/qc_summary.ipynb",
    conda:
        "../envs/model.yml"
    notebook:
        "notebooks/qc_summary.py.ipynb"
