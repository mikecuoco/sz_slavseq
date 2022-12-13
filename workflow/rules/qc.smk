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
    output:
        "{outdir}/results/multiqc.html",
    log:
        "{outdir}/results/multiqc.log",
    wrapper:
        "v1.21.0/bio/multiqc"
