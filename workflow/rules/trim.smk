rule cutadapt1:
    input:
        get_cutadapt_input,
    output:
        fastq1="{outdir}/results/cutadapt1/{donor}/{dna_type}/{sample}_R1.fastq.gz",
        fastq2="{outdir}/results/cutadapt1/{donor}/{dna_type}/{sample}_R2.fastq.gz",
        qc="{outdir}/results/cutadapt1/{donor}/{dna_type}/{sample}_qc.txt",
    params:
        adapters=f"--front={config['adapters']['r1_adapter']} --adapter={Seq(config['adapters']['nested_primer']).reverse_complement()}{Seq(config['adapters']['r2_adapter']).reverse_complement()}",
        extra="--minimum-length=36 --quality-base=33 --quality-cutoff=28 --overlap=5 --times=4",
    log:
        "{outdir}/results/cutadapt1/{donor}/{dna_type}/{sample}.log",
    threads: 2
    wrapper:
        "v1.19.2/bio/cutadapt/pe"


rule cutadapt2:
    input:
        [rules.cutadapt1.output.fastq1, rules.cutadapt1.output.fastq2],
    output:
        fastq1="{outdir}/results/cutadapt2/{donor}/{dna_type}/{sample}_R1.fastq.gz",
        fastq2="{outdir}/results/cutadapt2/{donor}/{dna_type}/{sample}_R2.fastq.gz",
        qc="{outdir}/results/cutadapt2/{donor}/{dna_type}/{sample}_qc.txt",
    params:
        adapters=f"--front={config['adapters']['r2_adapter']}{config['adapters']['nested_primer']} --adapter={Seq(config['adapters']['r1_adapter']).reverse_complement()}",
        extra="--minimum-length=36 --quality-base=33 --quality-cutoff=28 --overlap=5 --times=4",
    log:
        "{outdir}/results/cutadapt2/{donor}/{dna_type}/{sample}.log",
    threads: 2
    wrapper:
        "v1.19.2/bio/cutadapt/pe"
