rule bwa_index:
    input:
        rules.gen_ref.output[0],
    output:
        idx=multiext(
            f"{outdir}/resources/{{ref}}/{{ref}}{region_name}",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
       f"{outdir}/resources/{{ref}}/bwa_index.log",
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
        f"{outdir}/results/bwa_mem/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.bam",
    log:
        f"{outdir}/results/bwa_mem/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.log",
    params:
        extra="-T 19",
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sorting="samtools",
    threads: 4
    cache: True
    wrapper:
        "v1.18.0/bio/bwa/mem"


rule rmdup:
    input:
        rules.bwa_mem.output,
    output:
        f"{outdir}/results/rmdup/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.bam",
    log:
        f"{outdir}/results/rmdup/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.log",
    conda:
        "../envs/align.yml"
    script:
        "../scripts/slavseq_rmdup_hts.py"


rule install_gapafim:
    output:
        directory(f"{outdir}/resources/gapafim"),
    conda:
        "../envs/align.yml"
    log:
        f"{outdir}/resources/install_gapafim.log",
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        mkdir -p resources && cd resources
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
        f"{outdir}/results/tags/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.bam",
    log:
        f"{outdir}/results/tags/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.err",
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


rule tabix:
    input:
        rules.tags.output,
    output:
        bgz=f"{outdir}/results/tabix/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.bgz",
        tbi=f"{outdir}/results/tabix/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.bgz.tbi",
    log:
        f"{outdir}/results/tabix/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.log",
    conda:
        "../envs/align.yml"
    shell:
        """
        samtools view {input} | \
            workflow/scripts/sam_to_tabix.py | \
            sort --temporary-directory=results/tabix/{wildcards.sample} --buffer-size=10G -k1,1 -k2,2n -k3,3n | \
            bgzip -c > {output.bgz} 2> {log} 

        tabix -s 1 -b 2 -e 3 -0 {output.bgz} >> {log} 2>&1
        """
