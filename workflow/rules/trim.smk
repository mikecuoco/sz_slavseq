from Bio.Seq import Seq


def get_cutadapt_input(wildcards):
    sample = samples.loc[wildcards.donor, wildcards.sample]

    # for debugging
    # print("{} {} {}".format(wildcards.sample, sample["R1"], sample["R2"]))
    return [sample["R1"], sample["R2"]]


rule cutadapt:
    input:
        get_cutadapt_input,
    output:
        fastq1="{outdir}/results/trim/{donor}/{sample}_R1.fastq.gz",
        fastq2="{outdir}/results/trim/{donor}/{sample}_R2.fastq.gz",
        r1_qc="{outdir}/results/trim/{donor}/{sample}_R1.qc.txt",
        r2_qc="{outdir}/results/trim/{donor}/{sample}_R2.qc.txt",
    params:
        r1_front=config["adapters"]["r1_adapter"],
        r2_front=config["adapters"]["r2_adapter"] + config["adapters"]["nested_primer"],
        r1_end=Seq(
            config["adapters"]["r2_adapter"] + config["adapters"]["nested_primer"]
        ).reverse_complement(),
        r2_end=Seq(config["adapters"]["r1_adapter"]).reverse_complement(),
        extra="--minimum-length=36 --quality-base=33 --quality-cutoff=28 --overlap=5 --times=4",
    log:
        "{outdir}/results/trim/{donor}/{sample}.log",
    threads: 8
    conda:
        "../envs/trim.yml"
    shell:
        """
        touch {log} && exec > {log} 2>&1

        tmpdir=$(mktemp -d)

        cutadapt -j {threads} {params.extra} \
            --front={params.r1_front} --adapter={params.r1_end} \
            --paired-output $tmpdir/tmp.2.fq.gz -o $tmpdir/tmp.1.fq.gz \
            {input[0]} {input[1]} > {output.r1_qc}
        cutadapt -j {threads} {params.extra} \
            --front={params.r2_front} --adapter={params.r2_end} \
            --paired-output {output.fastq1} -o {output.fastq2} \
            $tmpdir/tmp.2.fq.gz $tmpdir/tmp.1.fq.gz > {output.r2_qc}

        rm -rf $tmpdir
        """
