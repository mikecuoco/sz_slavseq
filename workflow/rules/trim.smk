from Bio.Seq import Seq


# remove adapters and low quality reads
rule trim_adapters:
    input:
        lambda wc: [
            samples.loc[wc.donor, wc.sample]["R1"],
            samples.loc[wc.donor, wc.sample]["R2"],
        ],
    output:
        fastq1="{outdir}/results/fastq/{donor}/{sample}_R1.trimmed.fastq.gz",
        fastq2="{outdir}/results/fastq/{donor}/{sample}_R2.trimmed.fastq.gz",
        qc="{outdir}/results/fastq/{donor}/{sample}.trimmed.qc.txt",
    params:
        adapters="-a "
        + config["adapters"]["r1_adapter"]
        + " -A "
        + config["adapters"]["r2_adapter"]
        + " -g "
        + Seq(config["adapters"]["r2_adapter"]).reverse_complement()
        + " -G "
        + Seq(config["adapters"]["r1_adapter"]).reverse_complement(),
        extra="--minimum-length 36 --overlap 5 --times 4",
    log:
        "{outdir}/results/fastqc/{donor}/{sample}.trim_adapters.log",
    threads: 4
    wrapper:
        "v1.31.1/bio/cutadapt/pe"


# filter out read pairs where read2 doesn't contain the nested primer (that binds to LINE1)
rule filter_read2:
    input:
        [rules.trim_adapters.output.fastq1, rules.trim_adapters.output.fastq2],
    output:
        fastq1=rules.trim_adapters.output.fastq1.replace("trimmed", "filtered"),
        fastq2=rules.trim_adapters.output.fastq2.replace("trimmed", "filtered"),
        qc=rules.trim_adapters.output.qc.replace("trimmed", "filtered"),
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-G " + config["adapters"]["L1oligo_downstream"],
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--action=none --discard-untrimmed --pair-filter both --minimum-length 36 --quality-cutoff 20 --overlap 15",
    log:
        rules.trim_adapters.log[0].replace("trim_adapters", "filter_read2"),
    threads: 4
    wrapper:
        "v1.31.1/bio/cutadapt/pe"
