rule bwa_mem:
    input:
        idx=config["genome"]["bwa"] + ".amb",
        r1=rules.filter_l1hs.output.r1,
        r2=rules.filter_l1hs.output.r2,
        merged=rules.filter_l1hs.output.merged,
    output:
        "{outdir}/results/{genome}/align/{donor}/{sample}.genome.bam",
    log:
        bwa="{outdir}/results/{genome}/align/{donor}/{sample}.bwa.log",
        samblaster="{outdir}/results/{genome}/align/{donor}/{sample}.samblaster.log",
    threads: 2
    conda:
        "../envs/align.lock.yml"
    shell:
        """
        SAM=$(echo {output} | sed 's/.bam/.sam/g')
        trap "rm -f $SAM" EXIT

        # make sure r1 and r2 are the same length (for testing only)
        if [[ "{wildcards.genome}" == *"chr2122"* ]]; then
            NR1=$(zcat {input.r1} | wc -l)
            NR2=$(zcat {input.r2} | wc -l)
            if [ $NR1 -ne $NR2 ]; then
                echo "Error: r1 and r2 are not the same length"
                exit 1
            fi
        fi

        # run bwa mem and samblaster
        # do not add -T! will cause some mates to be lost
        IDX=$(echo {input.idx} | sed 's/.amb//g') # remove .amb from idx for bwa mem
        bwa mem -t {threads} -M $IDX {input.r1} {input.r2} 2> {log.bwa} | samblaster --addMateTags > $SAM 2> {log.samblaster}
        bwa mem -t {threads} -M $IDX {input.merged} 2>> {log.bwa} | samblaster --ignoreUnmated 2>> {log.samblaster} | samtools view >> $SAM
        samtools view -Sb -F 256 $SAM > {output}
        """


if not Path(config["genome"]["bwa"] + ".amb").exists():
    raise ValueError("BWA index not found for genome: " + config["genome"]["name"])


rule tag:
    input:
        genome_bam=rules.bwa_mem.output[0],
        line1_bam=rules.filter_l1hs.output.bam,
    output:
        rules.bwa_mem.output[0].replace("genome.", "tagged."),
    log:
        rules.bwa_mem.log.bwa.replace("bwa", "tag_reads"),
    conda:
        "../envs/align.lock.yml"
    script:
        "../scripts/tag_reads.py"


rule sort:
    input:
        rules.tag.output,
    output:
        rules.tag.output[0].replace("tagged", "tagged.sorted"),
    log:
        rules.tag.log[0].replace("tag_reads", "sort"),
    conda:
        "../envs/align.lock.yml"
    shell:
        "sambamba sort -o {output} {input} 2> {log}"


rule index:
    input:
        rules.sort.output,
    output:
        rules.sort.output[0].replace("bam", "bam.bai"),
    log:
        rules.sort.log[0].replace("sort", "index"),
    conda:
        "../envs/align.lock.yml"
    shell:
        "sambamba index {input} 2> {log}"


rule bigwig:
    input:
        rules.sort.output,
    output:
        r1=rules.sort.output[0].replace("bam", "r1.bw"),
        r2=rules.sort.output[0].replace("bam", "r2.bw"),
    log:
        rules.sort.log[0].replace("sort", "bigwig"),
    params:
        effective_genome_size=config["genome"]["size"],
        mapq=30,
    shell:
        """
        bamCoverage -b {input} -o {output.r1} -p {threads} \
            --minMappingQuality {params.mapq} \
            --effectiveGenomeSize {params.effective_genome_size} \
            --samFlagInclude 64
        bamCoverage -b {input} -o {output.r2} -p {threads} \
            --minMappingQuality {params.mapq} \
            --effectiveGenomeSize {params.effective_genome_size} \
            --samFlagInclude 128
        """


rule flagstat:
    input:
        rules.sort.output,
    output:
        rules.sort.output[0].replace("bam", "flagstat.txt"),
    log:
        rules.sort.log[0].replace("sort", "flagstat"),
    threads: 1
    conda:
        "../envs/align.lock.yml"
    shell:
        "sambamba flagstat {input} > {output} 2> {log}"


rule map:
    input:
        expand(
            expand(
                rules.flagstat.output,
                zip,
                sample=samples["sample_id"],
                donor=samples["donor_id"],
                allow_missing=True,
            ),
            genome=config["genome"]["name"],
            outdir=config["outdir"],
        ),
