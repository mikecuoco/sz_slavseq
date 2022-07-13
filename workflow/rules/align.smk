rule bwa_mem:
    input:
        reads=[rules.cutadapt.output.r1, rules.cutadapt.output.r2],
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
    conda: "../envs/perl.yml"
    shell:
        "workflow/scripts/slavseq_rmdup_hts.pl {input} {output}"