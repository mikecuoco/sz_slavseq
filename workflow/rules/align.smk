rule bwa_index:
    input:
        rules.gen_ref.output[0],
    output:
        idx=multiext(
            f"{{outdir}}/resources/{{ref}}/{{ref}}{region_name}",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "{outdir}/resources/{ref}/bwa_index.log",
    cache: True
    params:
        algorithm="bwtsw",
    wrapper:
        "v1.21.6/bio/bwa/index"


rule bwa_mem:
    input:
        reads=[rules.cutadapt.output.fastq1, rules.cutadapt.output.fastq2],
        idx=rules.bwa_index.output.idx,
    output:
        bam="{outdir}/results/align/bwa_mem/{ref}/{donor}/{dna_type}/{sample}.bam",
        index="{outdir}/results/align/bwa_mem/{ref}/{donor}/{dna_type}/{sample}.bam.bai",
    log:
        "{outdir}/results/align/bwa_mem/{ref}/{donor}/{dna_type}/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{donor}\tSM:{sample}'",
        samblaster_extra="-r --addMateTags",
    threads: 8
    wrapper:
        "v1.21.6/bio/bwa/mem-samblaster"


rule install_gapafim:
    output:
        directory("{outdir}/resources/gapafim"),
    conda:
        "../envs/tags.yml"
    log:
        "{outdir}/resources/install_gapafim.log",
    shell:
        """
        touch {log} && exec 1>{log} 2>&1
        mkdir -p $(dirname {output}) && cd $(dirname {output})
        git clone https://github.com/apuapaquola/gapafim.git
        cd gapafim/Gapafim
        perl Makefile.PL
        make
        make install
        """


rule tags:
    input:
        bam=rules.bwa_mem.output.bam,
        fa=rules.gen_ref.output[0],
        gapafim=rules.install_gapafim.output,
    output:
        bam="{outdir}/results/align/tags/{ref}/{donor}/{dna_type}/{sample}.bam",
        index="{outdir}/results/align/tags/{ref}/{donor}/{dna_type}/{sample}.bam.bai",
    log:
        "{outdir}/results/align/tags/{ref}/{donor}/{dna_type}/{sample}.err",
    conda:
        "../envs/tags.yml"
    params:
        consensus="ATGTACCCTAAAACTTAGAGTATAATAAA",
        prefix_length=len("ATGTACCCTAAAACTTAGAGTATAATAAA") + 2,
        r1_flank_length=750,
        r2_flank_length=len("ATGTACCCTAAAACTTAGAGTATAATAAA") + 2,
        soft_clip_length_threshold=5,
    shell:
        """
        touch {log} && exec 2> {log}
        (samtools sort -n {input.bam} | \
            samtools view -h | \
            workflow/scripts/add_tags_hts.pl \
                --genome_fasta_file {input.fa} \
                --prefix_length {params.prefix_length} \
                --consensus {params.consensus} \
                --r1_flank_length {params.r1_flank_length} \
                --r2_flank_length {params.r2_flank_length} \
                --soft_clip_length_threshold {params.soft_clip_length_threshold} | \
                samtools view -S -b - | \
                samtools sort - -o {output.bam})
        samtools index {output.bam}
        """
