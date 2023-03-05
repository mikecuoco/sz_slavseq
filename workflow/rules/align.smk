rule bwa_index:
    input:
        bwakit=rules.install_bwakit.output,
        fa=rules.gen_ref.output[0],
    output:
        idx=multiext(
            f"{{outdir}}/resources/hs38DH{region_name}.fa",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "{outdir}/resources/bwa_index.log",
    shell:
        "{input.bwakit}/bwa index {input.fa} > {log} 2>&1"


rule bwa_mem:
    input:
        bwakit=rules.install_bwakit.output,
        idx=rules.bwa_index.output.idx,
        fa=rules.gen_ref.output[0],
        reads=[rules.cutadapt.output.fastq1, rules.cutadapt.output.fastq2],
    output:
        "{outdir}/results/align/bwa_mem/{donor}/{dna_type}/{sample}.aln.bam",
    threads: 4
    shell:
        """
        prefix="$(dirname {output})/$(basename {output} .aln.bam)"
        idxbase="$(dirname {input.idx[0]})/$(basename {input.idx[0]} .amb)"

        # -s sort option doesn't work
        {input.bwakit}/run-bwamem \
            -R "@RG\\tID:{wildcards.donor}\\tSM:{wildcards.sample}\\tPL:ILLUMINA" \
            -d \
            -t {threads} \
            -o $prefix \
            $idxbase \
            {input.reads} | bash
        """


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
        bam=rules.bwa_mem.output,
        fa=rules.gen_ref.output[0],
        gapafim=rules.install_gapafim.output,
    output:
        "{outdir}/results/align/tags/{donor}/{dna_type}/{sample}.bam",
    log:
        "{outdir}/results/align/tags/{donor}/{dna_type}/{sample}.err",
    conda:
        "../envs/align.yml"
    params:
        consensus="ATGTACCCTAAAACTTAGAGTATAATAAA",
        prefix_length=len("ATGTACCCTAAAACTTAGAGTATAATAAA") + 2,
        r1_flank_length=750,
        r2_flank_length=len("ATGTACCCTAAAACTTAGAGTATAATAAA") + 2,
        soft_clip_length_threshold=5,
    shell:
        """
        touch {log} && exec 2>{log}

        (samtools sort -n {input.bam} | \
            samtools view -h | \
            workflow/scripts/add_tags_hts.pl \
                --genome_fasta_file {input.fa} \
                --prefix_length {params.prefix_length} \
                --consensus {params.consensus} \
                --r1_flank_length {params.r1_flank_length} \
                --r2_flank_length {params.r2_flank_length} \
                --soft_clip_length_threshold {params.soft_clip_length_threshold} | \
                samtools view -S -b - > {output})
        """


rule tabix:
    input:
        rules.tags.output,
    output:
        bgz="{outdir}/results/align/tabix/{donor}/{dna_type}/{sample}.bgz",
        tbi="{outdir}/results/align/tabix/{donor}/{dna_type}/{sample}.bgz.tbi",
    log:
        "{outdir}/results/align/tabix/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/align.yml"
    shell:
        """
        workflow/scripts/sam_to_tabix.py {input} | \
            sort --temporary-directory=results/tabix/{wildcards.sample} --buffer-size=10G -k1,1 -k2,2n -k3,3n | \
            bgzip -c > {output.bgz} 2> {log}

        tabix -s 1 -b 2 -e 3 -0 {output.bgz} >> {log} 2>&1
        """
