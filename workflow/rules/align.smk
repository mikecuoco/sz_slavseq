rule bwa_index:
    input:
        fa=f"{{outdir}}/resources/{genome_name}.fa",
    output:
        idx=multiext(
            f"{{outdir}}/resources/{genome_name}.fa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        f"{{outdir}}/resources/{genome_name}.bwa_index.log",
    conda:
        "../envs/align.yml"
    params:
        algorithm="bwtsw",
    shell:
        """
        bwa index -a {params.algorithm} {input.fa} > {log} 2>&1
        """


rule bwa_mem:
    input:
        idx=rules.bwa_index.output.idx,
        fa=rules.get_genome.output.fa,
        reads=[rules.filter_read2.output.fastq1, rules.filter_read2.output.fastq2],
    output:
        "{outdir}/results/align/{donor}/{sample}.bam",
    log:
        bwa="{outdir}/results/align/{donor}/{sample}.bwa.log",
        samblaster="{outdir}/results/align/{donor}/{sample}.samblaster_genome.log",
    threads: 4
    conda:
        "../envs/align.yml"
    params:
        min_as=19,  # minimum alignment score
    shell:
        """
        bwa mem -T {params.min_as} -t {threads} {input.fa} {input.reads} 2> {log.bwa} \
        | samblaster --addMateTags --removeDups 2> {log.samblaster} \
        | samtools view -Sb - > {output}
        """


rule sambamba_sort:
    input:
        rules.bwa_mem.output,
    output:
        rules.bwa_mem.output[0].replace(".bam", ".sorted.bam"),
    log:
        rules.bwa_mem.log.bwa.replace("bwa_genome", "sort"),
    threads: 4
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
