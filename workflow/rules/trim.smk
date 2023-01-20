from Bio.Seq import Seq


def get_cutadapt_input(wildcards):
    sample = samples.loc[wildcards.sample, wildcards.donor, wildcards.dna_type]

    if "R1" in sample:
        return [sample["R1"], sample["R2"]]
    else:
        accession = sample["sra"]
        return expand(
            "results/fastq/{accession}_{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )


rule cutadapt:
    input:
        get_cutadapt_input,
    output:
        fastq1="{outdir}/results/trim/{donor}/{dna_type}/{sample}_R1.fastq.gz",
        fastq2="{outdir}/results/trim/{donor}/{dna_type}/{sample}_R2.fastq.gz",
        r1_qc="{outdir}/results/trim/{donor}/{dna_type}/{sample}_R1.qc.txt",
        r2_qc="{outdir}/results/trim/{donor}/{dna_type}/{sample}_R2.qc.txt",
    params:
        r1_front=config["adapters"]["r1_adapter"],
        r2_front=config["adapters"]["r2_adapter"] + config["adapters"]["nested_primer"],
        r1_end=Seq(
            config["adapters"]["r2_adapter"] + config["adapters"]["nested_primer"]
        ).reverse_complement(),
        r2_end=Seq(config["adapters"]["r1_adapter"]).reverse_complement(),
        extra="--minimum-length=36 --quality-base=33 --quality-cutoff=28 --overlap=5 --times=4",
    log:
        "{outdir}/results/trim/{donor}/{dna_type}/{sample}.log",
    threads: 2
    conda:
        "../envs/trim.yml"
    shadow:
        "shallow"
    shell:
        """ 
        touch {log} && exec > {log} 2>&1

        cutadapt -j {threads} {params.extra} \
            --front={params.r1_front} --adapter={params.r1_end} \
            --paired-output tmp.2.{wildcards.sample}.fastq -o tmp.1.{wildcards.sample}.fastq \
            {input[0]} {input[1]} > {output.r1_qc}
        cutadapt -j {threads} {params.extra} \
            --front={params.r2_front} --adapter={params.r2_end} \
            --paired-output {output.fastq1} -o {output.fastq2} \
            tmp.2.{wildcards.sample}.fastq tmp.1.{wildcards.sample}.fastq > {output.r2_qc}

        rm -f tmp.2.{wildcards.sample}.fastq tmp.1.{wildcards.sample}.fastq
        """
