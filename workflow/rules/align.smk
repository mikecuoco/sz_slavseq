rule bwa_mem:
    input:
        reads=[rules.cutadapt2.output.fastq1, rules.cutadapt2.output.fastq2],
        idx=expand("resources/{ref}/genome.fa{ext}", ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], ref=config["ref"])
    output: "results/bwa_mem/{sample}/{donor}_{type}.bam"
    log: "results/bwa_mem/{sample}/{donor}_{type}.log"
    params: 
        extra="-T 19",
        sorting="samtools"
    threads: 4
    wrapper:
        "v1.7.0/bio/bwa/mem"

rule rmdup:
    input: rules.bwa_mem.output
    output: "results/rmdup/{sample}/{donor}_{type}.bam"
    log: "results/rmdup/{sample}/{donor}_{type}.log"
    conda: "../envs/env.yml"
    shell:
        "workflow/scripts/slavseq_rmdup_hts.pl {input} {output} > {log} 2>&1"

rule tags:
    input: 
        bam=rules.rmdup.output,
        ref=expand("resources/{ref}/genome.fa", ref=config["ref"])
    output: "results/tags/{sample}/{donor}_{type}.bam"
    log: "results/tags/{sample}/{donor}_{type}.err"
    conda: "../envs/env.yml"
    shell:
        '''
            # install gapafim
            # TODO: check if this will work when running in parallel
            if [ ! -d gapafim ]; then
                git clone https://github.com/apuapaquola/gapafim.git
                cd gapafim/Gapafim
                perl Makefile.PL
                make
                make install
                cd ../..
            fi

            # set inputs
            export CONSENSUS='ATGTACCCTAAAACTTAGAGTATAATAAA'
            PREFIX_LENGTH=`perl -e 'print length($ENV{{CONSENSUS}})+2'`
            R1_FLANK_LENGTH=750
            R2_FLANK_LENGTH=${{PREFIX_LENGTH}}
            SOFT_CLIP_LENGTH_THRESHOLD=5

            (samtools sort -n {input.bam} | \
                samtools view -h | \
                workflow/scripts/add_tags_hts.pl \
                    --genome_fasta_file {input.ref} \
                    --prefix_length ${{PREFIX_LENGTH}} \
                    --consensus ${{CONSENSUS}} \
                    --r1_flank_length ${{R1_FLANK_LENGTH}} \
                    --r2_flank_length ${{R2_FLANK_LENGTH}} \
                    --soft_clip_length_threshold ${{SOFT_CLIP_LENGTH_THRESHOLD}} | \
                    samtools view -S -b - > {output}) 2> {log}
        '''

rule tabix:
    input: rules.tags.output
    output: 
        bgz = "results/tabix/{sample}/{donor}_{type}.bgz",
        tbi = "results/tabix/{sample}/{donor}_{type}.bgz.tbi"
    log: "results/tabix/{sample}/{donor}_{type}.log"
    conda: "../envs/env.yml"
    shell:
        '''
        samtools view {input} | \
            workflow/scripts/sam_to_tabix.py | \
            sort --temporary-directory=results/tabix/{wildcards.sample} --buffer-size=10G -k1,1 -k2,2n -k3,3n | \
            bgzip -c > {output.bgz} 2> {log} 
        
        tabix -s 1 -b 2 -e 3 -0 {output.bgz} >> {log} 2>&1
        '''