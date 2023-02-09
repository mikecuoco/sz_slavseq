if config["genome"]["region"] != "all":
    fai = expand(
        rules.gen_ref.output[1], ref=config["genome"]["build"], outdir=config["outdir"]
    )[0]
    with open(fai, "r") as f:
        genome_size = sum([int(x.split("\t")[1]) for x in f.readlines()])

rule macs2:
    input:
        rules.filter.output
    output:
        peaks="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_peaks.narrowPeak",
        summits="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_summits.bed",
        xls="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_peaks.xls",
        cutoff="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_cutoff_analysis.txt",
    log:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/peaks.yml"
    params:
        genome_size="hs" if config["genome"]["region"] == "all" else genome_size,
        qValue_cutoff=0.05,
        extsize=200,
    shell:
        """
        macs2 callpeak \
            -g {params.genome_size} \
            -t {input} \
            -q {params.qValue_cutoff} \
            --nomodel \
            --extsize {params.extsize} \
            --format BAM \
            --name {wildcards.sample} \
            --outdir $(dirname {output.peaks}) \
            --cutoff-analysis 2> {log}
        """


rule peaks_report:
    input:
        bulk_peaks=lambda wc: expand(
            rules.macs2.output.peaks,
            sample=samples.loc[
                (samples["dna_type"] == "bulk")
                & (samples["donor_id"] == wc.donor),
            ]["sample_id"],
            dna_type="bulk",
            allow_missing=True,
        ),
        bulk_bam=lambda wc: expand(
            rules.tags.output.bam,
            sample=samples.loc[
                (samples["dna_type"] == "bulk")
                & (samples["donor_id"] == wc.donor),
            ]["sample_id"],
            dna_type="bulk",
            allow_missing=True,
        ),
        cell_peaks=lambda wc: expand(
            rules.macs2.output.peaks,
            sample=samples.loc[
                (samples["dna_type"] == "mda")
                & (samples["donor_id"] == wc.donor),
            ]["sample_id"],
            dna_type="mda",
            allow_missing=True,
        ),
        cell_bam=lambda wc: expand(
            rules.tags.output.bam,
            sample=samples.loc[
                (samples["dna_type"] == "mda")
                & (samples["donor_id"] == wc.donor),
            ]["sample_id"],
            dna_type="mda",
            allow_missing=True,
        ),
        rl1=rules.run_rmsk.output,
        knrgl=rules.get_donor_knrgl.output,
    output:
        "{outdir}/results/macs2/{ref}/{donor}/peak_report.ipynb",
    conda:
        "../envs/peaks.yml"
    log:
        notebook="{outdir}/results/macs2/{ref}/{donor}/peak_report.ipynb",
    notebook:
        "../notebooks/peaks_report.py.ipynb"


rule render_peaks_report:
    input:
        rules.peaks_report.output,
    output:
        "{outdir}/results/macs2/{ref}/{donor}/peak_report.html",
    conda:
        "../envs/peaks.yml"
    shell:
        "jupyter nbconvert --to html {input} --output $(basename {output})"
