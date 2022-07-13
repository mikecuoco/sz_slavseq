rule bwa_mem:
    input:
        reads=[rules.cutadapt2.output.r1, rules.cutadapt2.output.r2],
        idx=expand("resources/{ref}.fa{ext}", ext=[".amb", ".ann", ".bwt", ".pac", ".sa"], ref=config["ref"])
    output: "results/bwa_mam/{sample}/{donor}_{type}.bam"
    log: "results/bwa_mam/{sample}/{donor}_{type}.log"
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
        "workflow/scripts/slavseq_rmdup_hts.pl {input} {output}"

rule tags:
    input: 
        bam=rules.rmdup.output,
        ref=expand("resources/{ref}.fa", ref=config["ref"])
    output: "results/tags/{sample}/{donor}_{type}.bam"
    log: "results/tags/{sample}/{donor}_{type}.err"
    conda: "../envs/env.yml"
    shell:
        '''
            export CONSENSUS='ATGTACCCTAAAACTTAGAGTATAATAAA'
            PREFIX_LENGTH=`perl -e 'print length($ENV{{CONSENSUS}})+2'`
            R1_FLANK_LENGTH=750
            R2_FLANK_LENGTH=${{PREFIX_LENGTH}}
            SOFT_CLIP_LENGTH_THRESHOLD=5

            (samtools view -h {input.bam} | \
                workflow/scripts/add_tags_hts.pl \
                    --genome_fasta_file {input.ref} \
                    --prefix_length ${{PREFIX_LENGTH}} \
                    --consensus ${{CONSENSUS}} \
                    --r1_flank_length ${{R1_FLANK_LENGTH}} \
                    --r2_flank_length ${{R2_FLANK_LENGTH}} \
                    --soft_clip_length_threshold ${{SOFT_CLIP_LENGTH_THRESHOLD}} | \
                    samtools view -S -b - > {output}) 2> {log}
        '''