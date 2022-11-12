rule get_sra:
    output:
        f"{outdir}/results/fastq/{{accession}}_1.fastq.gz",
        f"{outdir}/results/fastq/{{accession}}_2.fastq.gz",
    log:
        f"{outdir}/results/fastq/{{accession}}.log",
    wrapper:
        "v1.10.0/bio/sra-tools/fasterq-dump"
