def get_bwa_input(wildcards):
    out = {}
    out["idx"] = config["genome"]["bwa"] + ".amb"
    if wildcards.reads == "trimmed":
        out["r1"] = rules.fastp.output.trimmed[0]
        out["r2"] = rules.fastp.output.trimmed[1]
        out["merged"] = rules.fastp.output.merged
    elif wildcards.reads == "filtered":
        out["r1"] = rules.filter_l1hs.output.r1
        out["r2"] = rules.filter_l1hs.output.r2
        out["merged"] = rules.filter_l1hs.output.merged
    else:
        raise ValueError("Invalid reads wildcard: " + wildcards.reads)
    return out


rule bwa_mem:
    input:
        unpack(get_bwa_input),
    output:
        "{outdir}/results/{genome}/{reads}/align/{donor}/{sample}.genome.bam",
    log:
        bwa="{outdir}/results/{genome}/{reads}/align/{donor}/{sample}.bwa.log",
        samblaster="{outdir}/results/{genome}/{reads}/align/{donor}/{sample}.samblaster.log",
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


rule reads_report:
    input:
        fastqc=expand(
            expand(
                rules.fastqc.output.zip,
                zip,
                donor=samples["donor_id"],
                sample=samples["sample_id"],
                allow_missing=True,
            ),
            zip,
            read=["R1", "R2", "R1", "R2", "merged", "R1", "R2", "merged"],
            fastq=[
                "raw",
                "raw",
                "trimmed",
                "trimmed",
                "trimmed",
                "filtered",
                "filtered",
                "filtered",
            ],
            allow_missing=True,
        ),
        flagstat=expand(
            expand(
                rules.flagstat.output,
                reads="filtered",
                allow_missing=True,
            ),
            zip,
            donor=samples["donor_id"],
            sample=samples["sample_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/{genome}/reads_report.ipynb",
    log:
        notebook="{outdir}/results/{genome}/reads_report.ipynb",
    conda:
        "../envs/model.yml"
    notebook:
        "../scripts/reads_report.py.ipynb"


# TODO: compute features here
rule coverage:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        l1hs_rmsk="resources/{genome}/{genome}.fasta.rmsk.l1hs.bed",
        chrom_sizes=config["genome"]["genome"],
        megane=expand(rules.merge_bed.output, vcf="megane", allow_missing=True),
    output:
        meg=temp("{outdir}/results/{genome}/{reads}/align/{donor}/{sample}.meg.bed"),
        cov_temp=temp(
            "{outdir}/results/{genome}/{reads}/align/{donor}/{sample}.coverage.bed.tmp"
        ),
        cov="{outdir}/results/{genome}/{reads}/align/{donor}/{sample}.coverage.bed",
    log:
        "{outdir}/results/{genome}/{reads}/align/{donor}/{sample}.coverage.log",
    conda:
        "../envs/align.yml"
    shell:
        """
        exec &>> {log}

        bedtools slop -i {input.megane} -g {input.chrom_sizes} -b 100 > {output.meg}
        samtools view -F 1412 -b {input.bam} -q 30 | bedtools coverage -a {output.meg} -b stdin -counts > {output.cov_temp}

        LIBSIZE=$(samtools view -F 1412 -c {input.bam})
        awk -v libsize=$LIBSIZE '{{print $0"\t"$(NF)/libsize * 1e6}}' {output.cov_temp} > {output.cov}
        """


rule coverage_report:
    input:
        cov=expand(
            rules.coverage.output.cov,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/{genome}/{reads}/align/coverage.ipynb",
    log:
        notebook="{outdir}/results/{genome}/{reads}/align/coverage.ipynb",
    conda:
        "../envs/model.yml"
    notebook:
        "../scripts/coverage.py.ipynb"


rule collect_coverage:
    input:
        expand(
            rules.reads_report.output,
            genome=config["genome"]["name"],
            outdir=config["outdir"],
        ),
        expand(
            rules.coverage_report.output,
            reads="filtered",
            genome=config["genome"]["name"],
            outdir=config["outdir"],
        ),
