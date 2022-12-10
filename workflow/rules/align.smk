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
        "v1.18.0/bio/bwa/index"


rule bwa_mem:
    input:
        reads=[rules.cutadapt2.output.fastq1, rules.cutadapt2.output.fastq2],
        idx=rules.bwa_index.output.idx,
    output:
        "{outdir}/results/bwa_mem/{ref}/{donor}/{dna_type}/{sample}.bam",
    log:
        "{outdir}/results/bwa_mem/{ref}/{donor}/{dna_type}/{sample}.log",
    params:
        extra="-T 19",  # Don’t output alignment with score lower than 19.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sorting="samtools",
    threads: 4
    cache: True
    wrapper:
        "v1.19.2/bio/bwa/mem"


rule rmdup:
    input:
        rules.bwa_mem.output,
    output:
        "{outdir}/results/rmdup/{ref}/{donor}/{dna_type}/{sample}.bam",
    log:
        "{outdir}/results/rmdup/{ref}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/align.yml"
    script:
        "../scripts/slavseq_rmdup_hts.py"


rule install_gapafim:
    output:
        directory("{outdir}/resources/gapafim"),
    conda:
        "../envs/align.yml"
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
        bam=rules.rmdup.output,
        fa=rules.gen_ref.output[0],
        gapafim=rules.install_gapafim.output,
    output:
        "{outdir}/results/tags/{ref}/{donor}/{dna_type}/{sample}.bam",
    log:
        "{outdir}/results/tags/{ref}/{donor}/{dna_type}/{sample}.err",
    conda:
        "../envs/align.yml"
    shell:
        """
        touch {log} && exec 2>{log} 

        # set inputs
        export CONSENSUS='ATGTACCCTAAAACTTAGAGTATAATAAA'
        PREFIX_LENGTH=`perl -e 'print length($ENV{{CONSENSUS}})+2'`
        R1_FLANK_LENGTH=750
        R2_FLANK_LENGTH=${{PREFIX_LENGTH}}
        SOFT_CLIP_LENGTH_THRESHOLD=5

        (samtools sort -n {input.bam} | \
            samtools view -h | \
            workflow/scripts/add_tags_hts.pl \
                --genome_fasta_file {input.fa} \
                --prefix_length ${{PREFIX_LENGTH}} \
                --consensus ${{CONSENSUS}} \
                --r1_flank_length ${{R1_FLANK_LENGTH}} \
                --r2_flank_length ${{R2_FLANK_LENGTH}} \
                --soft_clip_length_threshold ${{SOFT_CLIP_LENGTH_THRESHOLD}} | \
                samtools view -S -b - > {output}) 
        """


rule sort:
    input:
        rules.tags.output,
    output:
        "{outdir}/results/sort/{ref}/{donor}/{dna_type}/{sample}.bam",
    log:
        "{outdir}/results/sort/{ref}/{donor}/{dna_type}/{sample}.log",
    wrapper:
        "v1.19.2/bio/samtools/sort"


rule index:
    input:
        rules.sort.output,
    output:
        "{outdir}/results/sort/{ref}/{donor}/{dna_type}/{sample}.bam.bai",
    log:
        "{outdir}/results/sort/{ref}/{donor}/{dna_type}/{sample}.log",
    wrapper:
        "v1.19.2/bio/samtools/index"
