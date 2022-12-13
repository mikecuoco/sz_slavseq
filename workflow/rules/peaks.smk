rule macs2:
    input:
        rules.tags.output,
    output:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_peaks.narrowPeak",
    log:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.log"
    conda:
        "../envs/peaks.yml"
    shell:
        "macs2 callpeak -t {input} --outdir $(dirname {output})"
