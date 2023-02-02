rule peaks:
    input:
        bam=rules.tags.output.bam,
        index=rules.tags.output.index,
    output:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}.bed",
    conda:
        "../envs/peaks.yml"
    log:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}.log",
    script:
        "../scripts/peaks.py"


rule peaks_report:
    input:
        peaks=rules.peaks.output,
        rl1=rules.run_rmsk.output,
        knrgl=rules.get_donor_knrgl.output,
    output:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}.ipynb",
    conda:
        "../envs/peaks.yml"
    log:
        notebook="{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}.ipynb",
    notebook:
        "../notebooks/peaks_report.py.ipynb"


rule render_peaks_report:
    input:
        rules.peaks_report.output,
    output:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}.html",
    conda:
        "../envs/peaks.yml"
    shell:
        "jupyter nbconvert --to html {input} --output $(basename {output})"
