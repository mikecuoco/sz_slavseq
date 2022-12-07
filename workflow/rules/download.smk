rule get_sra:
    output:
        "{outdir}/results/fastq/{accession}_1.fastq.gz",
        "{outdir}/results/fastq/{accession}_2.fastq.gz",
    log:
        "{outdir}/results/fastq/{accession}.log",
    wrapper:
        "v1.10.0/bio/sra-tools/fasterq-dump"
