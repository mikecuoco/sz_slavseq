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
        "../envs/align.yml"
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
        "../envs/align.yml"
    script:
        "../scripts/tag_reads.py"


rule sort:
    input:
        rules.tag.output,
    output:
        rules.tag.output[0].replace("tagged", "tagged.sorted"),
    log:
        rules.tag.log[0].replace("tag_reads", "sort"),
    threads: 1
    wrapper:
        "v1.31.1/bio/sambamba/sort"


rule index:
    input:
        rules.sort.output,
    output:
        rules.sort.output[0].replace("bam", "bam.bai"),
    params:
        extra="",  # optional parameters
    log:
        rules.sort.log[0].replace("sort", "index"),
    threads: 1
    wrapper:
        "v1.31.1/bio/sambamba/index"


rule flagstat:
    input:
        rules.sort.output,
    output:
        rules.sort.output[0].replace("bam", "flagstat.txt"),
    log:
        rules.sort.log[0].replace("sort", "flagstat"),
    params:
        extra="",  # optional parameters
    threads: 1
    wrapper:
        "v3.3.6/bio/sambamba/flagstat"
