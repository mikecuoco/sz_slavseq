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
        rules.bwa_mem.output,
    output:
        "{outdir}/results/qc/flagstat/{donor}/{sample}.flagstat",
    log:
        "{outdir}/results/qc/flagstat/{donor}/{sample}.flagstat.log",
    wrapper:
        "v1.21.0/bio/samtools/flagstat"


rule depth:
    input:
        bams=rules.sambamba_sort.output,
        bai=rules.sambamba_index.output,
    output:
        "{outdir}/results/qc/depth/{donor}/{sample}.depth.txt",
    log:
        "{outdir}/results/qc/depth/{donor}/{sample}.depth.log",
    params:
        extra=" ",  # optional additional parameters as string
    wrapper:
        "v1.28.0/bio/samtools/depth"


rule reads_multiqc:
    input:
        expand(
            expand(
                rules.fastqc.output.html,
                fastq=["raw", "trimmed", "filtered"],
                read=["R1", "R2"],
                allow_missing=True,
            ),
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
        expand(
            rules.trim_adapters.output,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
        expand(
            rules.filter_read2.output,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/qc/multiqc_reads.html",
    log:
        "{outdir}/results/qc/multiqc_reads.log",
    params:
        extra='--title "SLAV-seq fastQC and cutadapt QC" --no-data-dir',
    wrapper:
        "v1.21.0/bio/multiqc"


rule aln_multiqc:
    input:
        expand(
            rules.flagstat.output,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
        expand(
            rules.depth.output,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/qc/multiqc.html",
    log:
        "{outdir}/results/qc/multiqc.log",
    params:
        extra='--title "SLAV-seq alignment QC" --no-data-dir',
    wrapper:
        "v1.21.0/bio/multiqc"
