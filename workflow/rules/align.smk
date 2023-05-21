rule bwa_index:
    input:
        fa="{outdir}/resources/{ref}.fa",
    output:
        idx=multiext(
            "{outdir}/resources/{ref}.fa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "{outdir}/resources/{ref}.bwa_index.log",
    conda:
        "../envs/align.yml"
    params:
        algorithm=lambda wc: "is" if "LINE" in wc.ref else "bwtsw",
    shell:
        """
        bwa index -a {params.algorithm} {input.fa} > {log} 2>&1
        """


rule bwa_mem_genome:
    input:
        idx=expand(rules.bwa_index.output.idx, ref=genome_name, allow_missing=True),
        fa=rules.get_genome.output.fa,
        reads=[rules.filter_read2.output.fastq1, rules.filter_read2.output.fastq2],
    output:
        "{outdir}/results/align/{donor}/{sample}.genome.bam",
    log:
        bwa="{outdir}/results/align/{donor}/{sample}.bwa_genome.log",
        samblaster="{outdir}/results/align/{donor}/{sample}.samblaster_genome.log",
    threads: 4
    conda:
        "../envs/align.yml"
    params:
        min_as=19,  # minimum alignment score
    shell:
        """
        bwa mem -T {params.min_as} -t {threads} {input.fa} {input.reads} 2> {log.bwa} \
        | samblaster --addMateTags 2> {log.samblaster} \
        | samtools view -Sb - > {output}
        """


rule bwa_mem_line1:
    input:
        idx=expand(rules.bwa_index.output.idx, ref="LINE1_lib", allow_missing=True),
        fa=rules.make_dfam_lib.output,
        reads=rules.filter_read2.output.fastq2,
    output:
        rules.bwa_mem_genome.output[0].replace("genome", "line1"),
    log:
        rules.bwa_mem_genome.log.bwa.replace("genome", "line1"),
    threads: 4
    conda:
        "../envs/align.yml"
    params:
        min_as=19,  # minimum alignment score
    shell:
        """
        bwa mem -T {params.min_as} -t {threads} {input.fa} {input.reads} 2> {log} > {output}
        """


rule tag_reads:
    input:
        genome_bam=rules.bwa_mem_genome.output[0],
        line1_bam=rules.bwa_mem_line1.output[0],
    output:
        rules.bwa_mem_genome.output[0].replace("genome", "tagged"),
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
