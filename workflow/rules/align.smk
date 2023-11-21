rule bwa_index:
    input:
        fa=rules.get_l1hs_consensus.output,
    output:
        idx=multiext(
            rules.get_l1hs_consensus.output[0],
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "resources/LINE1/bwa_index.log",
    conda:
        "../envs/align.yml"
    params:
        algorithm="is",
    shell:
        """
        bwa index -a {params.algorithm} {input.fa} > {log} 2>&1
        """


rule bwa_mem_genome:
    input:
        idx=config["genome"]["bwa"] + ".amb",
        reads=[rules.filter_read2.output.fastq1, rules.filter_read2.output.fastq2],
    output:
        temp("{outdir}/results/{genome}/align/{donor}/{sample}.genome.bam"),
    log:
        bwa="{outdir}/results/{genome}/align/{donor}/{sample}.bwa_genome.log",
        samblaster="{outdir}/results/{genome}/align/{donor}/{sample}.samblaster_genome.log",
    threads: 2
    conda:
        "../envs/align.yml"
    params:
        min_as=30,  # minimum alignment score (default = 30)
    shell:
        """
        # remove .amb from idx
        idx=$(echo {input.idx} | sed 's/.amb//g')
        bwa mem -T {params.min_as} -t {threads} $idx {input.reads} 2> {log.bwa} \
        | samblaster --addMateTags 2> {log.samblaster} \
        | samtools view -Sb - > {output}
        """


rule bwa_mem_line1:
    input:
        idx=rules.bwa_index.output.idx[0],
        fa=rules.get_l1hs_consensus.output,
        reads=rules.filter_read2.output.fastq2,
    output:
        temp(rules.bwa_mem_genome.output[0].replace("genome.", "line1.")),
    log:
        rules.bwa_mem_genome.log.bwa.replace("genome.", "line1."),
    threads: 1
    conda:
        "../envs/align.yml"
    params:
        min_as=30,  # minimum alignment score (default = 30)
    shell:
        """
        bwa mem -T {params.min_as} -t {threads} {input.fa} {input.reads} 2> {log} > {output}
        """


if not Path(config["genome"]["bwa"] + ".amb").exists():
    raise ValueError("BWA index not found for genome: " + config["genome"]["name"])


rule flagstat:
    input:
        "{outdir}/results/{genome}/align/{donor}/{sample}.genome.bam",
    output:
        "{outdir}/results/{genome}/align/{donor}/{sample}.genome.flagstat",
    log:
        "{outdir}/results/{genome}/align/{donor}/{sample}.genome.flagstat.log",
    wrapper:
        "v1.21.0/bio/samtools/flagstat"


rule tag_reads:
    input:
        genome_bam=rules.bwa_mem_genome.output[0],
        line1_bam=rules.bwa_mem_line1.output[0],
    output:
        temp(rules.bwa_mem_genome.output[0].replace("genome.", "tagged.")),
    log:
        rules.bwa_mem_genome.log.bwa.replace("bwa_genome", "tag_reads"),
    conda:
        "../envs/align.yml"
    script:
        "../scripts/tag_reads.py"


rule sambamba_sort:
    input:
        rules.tag_reads.output,
    output:
        rules.tag_reads.output[0].replace("tagged", "tagged.sorted"),
    log:
        rules.bwa_mem_genome.log.bwa.replace("bwa_genome", "sort"),
    threads: 1
    wrapper:
        "v1.31.1/bio/sambamba/sort"


rule sambamba_index:
    input:
        rules.sambamba_sort.output,
    output:
        rules.sambamba_sort.output[0].replace("bam", "bam.bai"),
    params:
        extra="",  # optional parameters
    log:
        rules.sambamba_sort.log[0].replace("sort", "index"),
    threads: 1
    wrapper:
        "v1.31.1/bio/sambamba/index"
