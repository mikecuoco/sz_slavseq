def get_fastqc_input(wildcards):
    if wildcards.trim == "none":
        reads = get_cutadapt_input(wildcards)
    if wildcards.trim == "cutadapt":
        reads = rules.cutadapt.output
    return reads[0] if wildcards.read == "R1" else reads[1]


rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="{outdir}/results/qc/fastqc/{donor}/{sample}_{trim}_{read}.html",
        zip="{outdir}/results/qc/fastqc/{donor}/{sample}_{trim}_{read}_fastqc.zip",
    log:
        "{outdir}/results/qc/fastqc/{donor}/{sample}_{trim}_{read}.log",
    wrapper:
        "v1.21.0/bio/fastqc"


rule flagstat:
    input:
        "{outdir}/results/align/{donor}/{sample}.{ref}.bam",
    output:
        "{outdir}/results/qc/flagstat/{donor}/{sample}.{ref}.flagstat",
    log:
        "{outdir}/results/qc/flagstat/{donor}/{sample}.{ref}.flagstat.log",
    wrapper:
        "v1.21.0/bio/samtools/flagstat"


rule depth:
    input:
        bams=rules.sambamba_sort.output,
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
                trim=["none", "cutadapt"],
                read=["R1", "R2"],
                allow_missing=True,
            ),
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
        expand(
            rules.cutadapt.output,
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
        extra='--config config/multiqc_config.yml --title "SLAV-seq reads" --no-data-dir',
    wrapper:
        "v1.21.0/bio/multiqc"


rule aln_multiqc:
    input:
        expand(
            expand(
                rules.flagstat.output,
                ref=["genome", "line1"],
                allow_missing=True,
            ),
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
        extra='--config config/multiqc_config.yml --title "SLAV-seq" --no-data-dir',
    wrapper:
        "v1.21.0/bio/multiqc"
