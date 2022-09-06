rule get_sra:
    output:
        "results/fastq/{accession}_1.fastq.gz",
        "results/fastq/{accession}_2.fastq.gz",
    log:
        "results/fastq/{accession}.log",
    wrapper:
        "v1.10.0/bio/sra-tools/fasterq-dump"
