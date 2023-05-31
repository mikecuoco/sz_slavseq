from Bio.Seq import Seq


# Remove adapters from read through in short fragments
# If r1 is trimmed, it has read through LINE1 into the r2 adapter
# Use --minimum-length 80:20 to ensure r1 is still long enough to map uniquely
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
        + Seq(config["adapters"]["r2"]).reverse_complement()
        + " -A "
        + Seq(config["adapters"]["r1"]).reverse_complement(),
        extra="--minimum-length 80 --quality-cutoff 30",
    log:
        "{outdir}/results/fastq/{donor}/{sample}.trim_adapters.log",
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
        adapters="-G "
        + config["adapters"]["L1_nested"]
        + config["adapters"]["L1_downstream"],
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        # --action=none --discard-untrimmed, don't trim, just filter
        # --times 2, try to remove adapter twice
        extra="--action=none --discard-untrimmed --pair-filter both --overlap 50 --times 2",
    log:
        rules.trim_adapters.log[0].replace("trim_adapters", "filter_read2"),
    threads: 4
    wrapper:
        "v1.31.1/bio/cutadapt/pe"
