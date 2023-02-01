rule bwa_index:
    input:
        rules.gen_ref.output[0],
    output:
        idx=multiext(
            f"{{outdir}}/resources/{{ref}}/{{ref}}{region_name}",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "{outdir}/resources/{ref}/bwa_index.log",
    cache: True
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.6/bio/bwa/index"


rule bwa_mem:
    input:
        reads=[rules.cutadapt.output.fastq1, rules.cutadapt.output.fastq2],
        idx=rules.bwa_index.output.idx,
    output:
        bam="{outdir}/results/align/bwa_mem/{ref}/{donor}/{dna_type}/{sample}.bam",
        index="{outdir}/results/align/bwa_mem/{ref}/{donor}/{dna_type}/{sample}.bam.bai",
    log:
        "{outdir}/results/align/bwa_mem/{ref}/{donor}/{dna_type}/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{donor}\tSM:{sample}'",
        samblaster_extra="-r --addMateTags",
    threads: 32
    wrapper:
        "v1.21.6/bio/bwa/mem-samblaster"
