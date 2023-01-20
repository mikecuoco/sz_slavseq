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
        "macs2 callpeak -t {input} --name {wildcards.sample} --outdir $(dirname {output}) --nomodel -extsize 750"

rule macs2_evaluate:
    input:
        rules.macs2.output,
        rules.tags.output,
        rules.get_labels.output,
    output:
        "{outdir}/results/macs2_eval/{ref}_{db}/{donor}/{dna}/{sample}.ipynb"
    log:
        notebook="{outdir}/results/macs2_eval/{ref}_{db}/{donor}/{dna}/{sample}.ipynb"
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/evaluate_macs2.py.ipynb"