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
        html="{outdir}/results/qc/fastqc/{donor}/{dna_type}/{sample}_{trim}_{read}.html",
        zip="{outdir}/results/qc/fastqc/{donor}/{dna_type}/{sample}_{trim}_{read}_fastqc.zip",
    log:
        "{outdir}/results/qc/fastqc/{donor}/{dna_type}/{sample}_{trim}_{read}.log",
    wrapper:
        "v1.21.0/bio/fastqc"


rule flagstat:
    input:
        rules.bwa_mem.output,
    output:
        "{outdir}/results/qc/flagstat/bwa_mem/{donor}/{dna_type}/{sample}.flagstat",
    log:
        "{outdir}/results/qc/flagstat/bwa_mem/{donor}/{dna_type}/{sample}.flagstat.log",
    wrapper:
        "v1.21.0/bio/samtools/flagstat"


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
            sample=samples["sample"],
            donor=samples["donor"],
            dna_type=samples["dna_type"],
            allow_missing=True,
        ),
        expand(
            rules.cutadapt.output,
            zip,
            sample=samples["sample"],
            donor=samples["donor"],
            dna_type=samples["dna_type"],
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
                allow_missing=True,
            ),
            zip,
            sample=samples["sample"],
            donor=samples["donor"],
            dna_type=samples["dna_type"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/qc/multiqc.html",
    log:
        "{outdir}/results/qc/multiqc.log",
    params:
        extra=lambda wildcards: f'--config config/multiqc_config.yml --title "SLAV-seq" --no-data-dir',
    wrapper:
        "v1.21.0/bio/multiqc"
