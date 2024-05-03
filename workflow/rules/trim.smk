from Bio.Seq import Seq


# Remove adapters from read through in short fragments
# If r1 is trimmed, it has read through LINE1 into the r2 adapter
rule fastp:
    input:
        sample=lambda wc: [
            samples.loc[wc.donor, wc.sample]["R1"],
            samples.loc[wc.donor, wc.sample]["R2"],
        ],
    output:
        trimmed=[
            "{outdir}/results/fastp/{donor}/{sample}_R1.fastq.gz",
            "{outdir}/results/fastp/{donor}/{sample}_R2.fastq.gz",
        ],
        merged="{outdir}/results/fastp/{donor}/{sample}_merged.fastq.gz",
        html="{outdir}/results/fastp/{donor}/{sample}.html",
        json="{outdir}/results/fastp/{donor}/{sample}.json",
    log:
        "{outdir}/results/fastp/{donor}/{sample}.log",
    params:
        adapters="--adapter_sequence "
        + Seq(config["adapters"]["r2"]).reverse_complement()
        + " --adapter_sequence_r2 "
        + Seq(config["adapters"]["r1"]).reverse_complement(),
        extra="-F 1 --length_required 30 --qualified_quality_phred 20 --dont_eval_duplication --correction --merge",
    threads: 2
    conda:
        "../envs/align.lock.yml"
    shell:
        """
        exec &>> {log}

        merged=$(mktemp --suffix=.fastq)
        trap "rm -f $merged" EXIT
        fastp --in1 {input.sample[0]} --in2 {input.sample[1]} \
            --out1 {output.trimmed[0]} --out2 {output.trimmed[1]} --merged_out $merged \
            --html {output.html} --json {output.json} \
            --thread {threads} {params.adapters} {params.extra}

        # remove merged entries less than 100bp
        seqtk seq -L 100 $merged | gzip > {output.merged}
        """


rule index_line1_consensus:
    input:
        fa=rules.get_line1_consensus.output.fa,
    output:
        idx=multiext(
            rules.get_line1_consensus.output.fa,
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "resources/LINE1/bwa_index.log",
    conda:
        "../envs/align.lock.yml"
    params:
        algorithm="is",
    shell:
        """
        bwa index -a {params.algorithm} {input.fa} > {log} 2>&1
        """


# only keep read pairs that map the best to L1HS
rule filter_l1hs:
    input:
        r1=rules.fastp.output.trimmed[0],
        r2=rules.fastp.output.trimmed[1],
        merged=rules.fastp.output.merged,
        fa=rules.get_line1_consensus.output.fa,
        idx=rules.index_line1_consensus.output.idx,
    output:
        bam="{outdir}/results/l1hs_filter/{donor}/{sample}.bam",
        r1="{outdir}/results/l1hs_filter/{donor}/{sample}_R1.fastq.gz",
        r2="{outdir}/results/l1hs_filter/{donor}/{sample}_R2.fastq.gz",
        merged="{outdir}/results/l1hs_filter/{donor}/{sample}_merged.fastq.gz",
    log:
        "{outdir}/results/l1hs_filter/{donor}/{sample}.log",
    conda:
        "../envs/align.lock.yml"
    threads: 2
    params:
        min_mapq=12,
    shell:
        """
        exec &>> {log}

        # count reads
        n_reads_in=$(zcat {input.r2} | wc -l)
        n_reads_in=$((n_reads_in + $(zcat {input.merged} | wc -l)))
        n_reads_in=$((n_reads_in / 4))

        # create temp files
        reads=$(mktemp)
        trap "rm -f $reads" EXIT

        # map + filter read2, only keep reads that map to L1HS with MAPQ >= params.min_mapq
        IDX=$(echo {input.idx[0]} | sed 's/.amb//g')
        # bwa mem -t {threads} $IDX <(zcat {input.r2} {input.merged}) | samtools view -h -F 256 > {output.bam}
        bwa mem -t {threads} $IDX <(zcat {input.r2} {input.merged}) | \
            awk '{{if ($1 ~ /^@/ || $3 == "L1HS_3end") print $0}}' | \
            samtools view -h -F 256 -q {params.min_mapq} -bS - > {output.bam}
        samtools view {output.bam} | awk '{{print $1}}' > $reads
        zcat {input.r1} | grep -A 3 -F -f $reads --no-group-separator - | gzip > {output.r1}
        zcat {input.r2} | grep -A 3 -F -f $reads --no-group-separator - | gzip > {output.r2}
        zcat {input.merged} | grep -A 3 -F -f $reads --no-group-separator - | gzip > {output.merged}

        # count reads
        n_reads_out=$(zcat {output.r2} | wc -l)
        n_reads_out=$((n_reads_out + $(zcat {output.merged} | wc -l)))
        n_reads_out=$((n_reads_out / 4))

        # write stats
        echo "Reads retained: $n_reads_out/$n_reads_in"
        """


def get_fastqc_input(wc):
    if wc.fastq == "raw":
        reads = [
            samples.loc[wc.donor, wc.sample][["R1"]],
            samples.loc[wc.donor, wc.sample][["R2"]],
        ]
    elif wc.fastq == "trimmed":
        reads = [rules.fastp.output.fastq1, rules.fastp.output.fastq2]
    elif wc.fastq == "filtered":
        reads = [rules.filter_l1hs.output.fastq1, rules.filter_l1hs.output.fastq2]
    return reads[0] if wc.read == "R1" else reads[1]


rule fastqc:
    input:
        get_fastqc_input,
    output:
        html="{outdir}/results/fastqc/{donor}/{sample}_{read}.{fastq}.html",
        zip="{outdir}/results/fastqc/{donor}/{sample}_{read}.{fastq}.fastqc.zip",
    log:
        "{outdir}/results/fastqc/{donor}/{sample}_{read}.{fastq}.log",
    wrapper:
        "v1.21.0/bio/fastqc"
