rule get_sra:
    output:
        "results/fastq/{accession}_1.fastq",
        "results/fastq/{accession}_2.fastq",
    log:
        "results/fastq/{accession}.log",
    wrapper:
        "0.77.0/bio/sra-tools/fasterq-dump"

rule cutadapt1:
    input:
        get_cutadapt_input
    output:
        fastq1="results/cutadapt1/{donor}/{dna_type}/{sample}_R1.fastq.gz",
        fastq2="results/cutadapt1/{donor}/{dna_type}/{sample}_R2.fastq.gz",
        qc="results/cutadapt1/{donor}/{dna_type}/{sample}_qc.txt",
    params:
        adapters=f"--front={config['adapters']['r1_adapter']} --adapter={Seq(config['adapters']['nested_primer']).reverse_complement()}{Seq(config['adapters']['r2_adapter']).reverse_complement()}",
        extra="--minimum-length=36 --quality-base=33 --quality-cutoff=28 --overlap=5 --times=4",
    log:
        "results/cutadapt1/{donor}/{dna_type}/{sample}.log",
    threads: 4
    conda:
        "../envs/env.yml"
    wrapper:
        "v1.7.1/bio/cutadapt/pe"


rule cutadapt2:
    input:
        [rules.cutadapt1.output.fastq1, rules.cutadapt1.output.fastq2],
    output:
        fastq1="results/cutadapt2/{donor}/{dna_type}/{sample}_R1.fastq.gz",
        fastq2="results/cutadapt2/{donor}/{dna_type}/{sample}_R2.fastq.gz",
        qc="results/cutadapt2/{donor}/{dna_type}/{sample}_qc.txt",
    params:
        adapters=f"--front={config['adapters']['r2_adapter']}{config['adapters']['nested_primer']} --adapter={Seq(config['adapters']['r1_adapter']).reverse_complement()}",
        extra="--minimum-length=36 --quality-base=33 --quality-cutoff=28 --overlap=5 --times=4",
    log:
        "results/cutadapt2/{donor}/{dna_type}/{sample}.log",
    threads: 4
    conda:
        "../envs/env.yml"
    wrapper:
        "v1.7.1/bio/cutadapt/pe"
