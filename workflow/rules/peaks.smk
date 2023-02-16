rule peaks:
    input:
        bam=rules.tags.output.bam,
        index=rules.tags.output.index,
    output:
        peaks="{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}.bed",
    conda:
        "../envs/peaks.yml"
    log:
        "{outdir}/results/peaks/{ref}/{donor}/{dna_type}/{sample}.log",
    script:
        "../scripts/peaks.py"


rule peaks_report:
    input:
        bulk_peaks=lambda wc: expand(
            rules.peaks.output.peaks,
            sample=samples.loc[
                (samples["dna_type"] == "bulk") & (samples["donor_id"] == wc.donor),
            ]["sample_id"],
            dna_type="bulk",
            allow_missing=True,
        ),
        bulk_bam=lambda wc: expand(
            rules.tags.output.bam,
            sample=samples.loc[
                (samples["dna_type"] == "bulk") & (samples["donor_id"] == wc.donor),
            ]["sample_id"],
            dna_type="bulk",
            allow_missing=True,
        ),
        cell_peaks=lambda wc: expand(
            rules.peaks.output.peaks,
            sample=samples.loc[
                (samples["dna_type"] == "mda") & (samples["donor_id"] == wc.donor),
            ]["sample_id"],
            dna_type="mda",
            allow_missing=True,
        ),
        cell_bam=lambda wc: expand(
            rules.tags.output.bam,
            sample=samples.loc[
                (samples["dna_type"] == "mda") & (samples["donor_id"] == wc.donor),
            ]["sample_id"],
            dna_type="mda",
            allow_missing=True,
        ),
        rl1=rules.run_rmsk.output,
        knrgl=rules.get_donor_knrgl.output,
    output:
        notebook="{outdir}/results/peaks/{ref}/{donor}/peaks_report.ipynb",
        bulk_insertions="{outdir}/results/peaks/{ref}/{donor}/bulk_insertions.tsv",
    conda:
        "../envs/peaks.yml"
    log:
        notebook="{outdir}/results/peaks/{ref}/{donor}/peaks_report.ipynb",
    notebook:
        "../notebooks/peaks_report.py.ipynb"


rule render_peaks_report:
    input:
        rules.peaks_report.output.notebook,
    output:
        "{outdir}/results/peaks/{ref}/{donor}/peak_report.html",
    conda:
        "../envs/peaks.yml"
    shell:
        "jupyter nbconvert --to html {input} --output $(basename {output})"
