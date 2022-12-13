def get_fastqc_input(wildcards):
    if wildcards.trim == "none":
        reads = get_cutadapt_input(wildcards)
    if wildcards.trim == "cutadapt1":
        reads = rules.cutadapt1.output
    if wildcards.trim == "cutadapt2":
        reads = rules.cutadapt2.output
    return reads[0] if wildcards.read == "R1" else reads[1]


rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="{outdir}/results/fastqc/{donor}/{dna_type}/{sample}_{trim}_{read}.html",
        zip="{outdir}/results/fastqc/{donor}/{dna_type}/{sample}_{trim}_{read}_fastqc.zip",
    log:
        "{outdir}/fastqc/{donor}/{dna_type}/{sample}_{trim}_{read}.log",
    wrapper:
        "v1.21.0/bio/fastqc"


rule flagstat:
    input:
        "{outdir}/results/{stage}/{ref}/{donor}/{dna_type}/{sample}.bam",
    output:
        "{outdir}/results/{stage}/{ref}/{donor}/{dna_type}/{sample}.flagstat",
    log:
        "{outdir}/results/{stage}/{ref}/{donor}/{dna_type}/{sample}.flagstat.log",
    wrapper:
        "v1.21.0/bio/samtools/flagstat"

# WIP
# rule feature_label_metrics:
#     input:
#         expand(rules.get_labels.output, donor=set(samples["donor"])),
#     output:
#         res="{outdir}/results/model/class_abundance.tsv",
#         plot="{outdir}/results/model/class_abundance.pdf",
#     log:
#         "{outdir}/results/class_abundance.log",
#     conda:
#         "../envs/model.yml"
#     script:
#         "../scripts/feature_label_metrics.py"


rule multiqc:
    input:
        expand(
            expand(
                rules.fastqc.output.html,
                trim=["none", "cutadapt1", "cutadapt2"],
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
            rules.cutadapt1.output,
            zip,
            sample=samples["sample"],
            donor=samples["donor"],
            dna_type=samples["dna_type"],
            allow_missing=True,
        ),
        expand(
            rules.cutadapt2.output,
            zip,
            sample=samples["sample"],
            donor=samples["donor"],
            dna_type=samples["dna_type"],
            allow_missing=True,
        ),
        expand(
            expand(
                rules.flagstat.output,
                stage=["bwa_mem", "rmdup"],
                allow_missing=True,
            ),
            zip,
            sample=samples["sample"],
            donor=samples["donor"],
            dna_type=samples["dna_type"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/{ref}_multiqc.html",
    log:
        "{outdir}/results/{ref}_multiqc.log",
    params:
        extra=lambda wildcards: f"--config config/multiqc_config.yml --title \"SLAV-seq  {wildcards.ref}\" --no-data-dir"
    wrapper:
        "v1.21.0/bio/multiqc"
