rule peaks:
    input:
        bam=rules.bwa_mem.output.bam,
        bai=rules.bwa_mem.output.index,
    output:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}_{t}.bed",
    conda:
        "../envs/peaks.yml"
    params:
        cluster_threshold=lambda wc: int(wc.t),  # must be >=1 or None
    log:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}_{t}.log",
    script:
        "../scripts/peaks.py"


rule peaks_report:
    input:
        peaks=rules.peaks.output,
        rl1=rules.run_rmsk.output,
        knrgl=rules.get_donor_knrgl.output,
    output:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}_cluster{t}.ipynb",
    conda:
        "../envs/peaks.yml"
    log:
        notebook="{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}_cluster{t}.ipynb",
    notebook:
        "../notebooks/peaks_report.py.ipynb"


rule compare_peaks_report:
    input:
        peaks=expand(
            rules.peaks.output,
            t=[0, 1, 2, 4, 8, 16, 32, 64, 128, 256],
            allow_missing=True,
        ),
        rl1=rules.run_rmsk.output,
        knrgl=rules.get_donor_knrgl.output,
    output:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}_compare.ipynb",
    conda:
        "../envs/peaks.yml"
    params:
        t=[0, 1, 2, 4, 8, 16, 32, 64, 128, 256],
    log:
        notebook="{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}_compare.ipynb",
    notebook:
        "../notebooks/compare_peaks_report.py.ipynb"


rule render_peaks_report:
    input:
        rules.peaks_report.output,
    output:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}_cluster{t}.html",
    conda:
        "../envs/peaks.yml"
    shell:
        "jupyter nbconvert --to html {input} --output {output}"
