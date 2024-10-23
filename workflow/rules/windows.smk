rule make_windows:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        chrom_sizes=config["genome"]["genome"],
    output:
        bed="{outdir}/results/{genome}/{reads}/windows/{donor}/{sample}.bed.gz",
        tabix="{outdir}/results/{genome}/{reads}/windows/{donor}/{sample}.bed.gz.tbi",
    log:
        "{outdir}/results/{genome}/{reads}/windows/{donor}/{sample}.make_windows.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/make_windows.py"


label_out = {}
for label in [
    "bulk",
    "xtea",
    "megane_gaussian",
    "megane_percentile",
    "megane_breakpoints",
    "graffite",
    "primer_sites",
    "l1hs",
    "l1pa2",
    "l1pa3",
    "l1pa4",
    "l1pa5",
    "l1pa6",
    "polyA",
    "polyT",
]:
    label_out[label] = rules.make_windows.output.bed.replace(
        ".bed.gz", f".{label}.bed.gz"
    )


# TODO: add closest and overlap labels
rule label_windows:
    input:
        **rmsk_anno,
        **knrgl_anno,
        bulk_regions=lambda wc: expand(
            rules.make_peaks.output.bed,
            sample=bulk.loc[wc.donor, "sample_id"],
            allow_missing=True,
        ),
        cell_regions=rules.make_windows.output.bed,
        primer_sites=rules.blast_primers.output.bed,
        fai=config["genome"]["fai"],
    output:
        **label_out,
    conda:
        "../envs/ref.lock.yml"
    log:
        rules.make_windows.log[0].replace("make_windows", "label"),
    shell:
        """
        exec &>> {log}

        bedtools intersect -a {input.bulk_regions} -b {input.cell_regions} -wao | bgzip -ci -I {output.bulk}.gzi > {output.bulk}
        bedtools intersect -a {input.primer_sites} -b {input.cell_regions} -wao | bgzip -ci -I {output.primer_sites}.gzi > {output.primer_sites}
        bedtools intersect -a {input.xtea} -b {input.cell_regions} -wao | bgzip -ci -I {output.xtea}.gzi > {output.xtea}
        bedtools intersect -a {input.megane_gaussian}  -b {input.cell_regions} -wao | bgzip -ci -I {output.megane_gaussian}.gzi > {output.megane_gaussian}
        bedtools intersect -a {input.megane_percentile} -b {input.cell_regions} -wao | bgzip -ci -I {output.megane_percentile}.gzi > {output.megane_percentile}
        bedtools intersect -a {input.megane_breakpoints} -b {input.cell_regions} -wao | bgzip -ci -I {output.megane_breakpoints}.gzi > {output.megane_breakpoints}
        bedtools intersect -a {input.graffite} -b {input.cell_regions} -wao | bgzip -ci -I {output.graffite}.gzi > {output.graffite}
        bedtools intersect -a {input.l1hs} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1hs}.gzi > {output.l1hs}
        bedtools intersect -a {input.l1pa2} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa2}.gzi > {output.l1pa2}
        bedtools intersect -a {input.l1pa3} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa3}.gzi > {output.l1pa3}
        bedtools intersect -a {input.l1pa4} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa4}.gzi > {output.l1pa4}
        bedtools intersect -a {input.l1pa5} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa5}.gzi > {output.l1pa5}
        bedtools intersect -a {input.l1pa6} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa6}.gzi > {output.l1pa6}
        bedtools intersect -a {input.polyA} -b {input.cell_regions} -wao | bgzip -ci -I {output.polyA}.gzi > {output.polyA}
        bedtools intersect -a {input.polyT} -b {input.cell_regions} -wao | bgzip -ci -I {output.polyT}.gzi > {output.polyT}
        """


rule merge_window_labels:
    input:
        regions=rules.make_windows.output.bed,
        annotations=expand(
            rules.label_windows.output.bulk.replace("bulk", "{annotation}"),
            annotation=label_out,
            allow_missing=True,
        ),
    output:
        rules.label_windows.output.bulk.replace("bulk.bed.gz", "labelled.bed.gz"),
    log:
        rules.label_windows.log[0].replace("label", "merge"),
    conda:
        "../envs/features.yml"
    script:
        "../scripts/merge_labels.py"


rule windows_report:
    input:
        cells=expand(
            rules.merge_window_labels.output,
            zip,
            donor=cells["donor_id"],
            sample=cells["sample_id"],
            allow_missing=True,
        ),
        bulk=rules.bulk_peaks_report.output.bed,
        meta=config["donors"],
    output:
        "{outdir}/results/{genome}/{reads}/windows/data.bed.gz",
    conda:
        "../envs/model.yml"
    log:
        notebook="{outdir}/results/{genome}/{reads}/windows/report.ipynb",
    notebook:
        "../scripts/windows_report.py.ipynb"


rule model_windows:
    input:
        data=rules.windows_report.output,
    output:
        model=rules.windows_report.output[0].replace("data.bed.gz", "model.pkl"),
        best_hp=rules.windows_report.output[0].replace("data.bed.gz", "best_hp.json"),
        history=rules.windows_report.output[0].replace("data.bed.gz", "history.log"),
        predictions=rules.windows_report.output[0].replace(
            "data.bed.gz", "predictions.pqt"
        ),
        shap_values=rules.windows_report.output[0].replace(
            "data.bed.gz", "shap_values.pqt"
        ),
    conda:
        "../envs/model_gpu.yml" if shutil.which("nvidia-smi") else "../envs/model.yml"
    params:
        max_iter=50,
        random_state=1,
        features=[
            "n_reads",
            "n_duplicates",
            "n_unique_5end",
            "three_end_clippedA_mean",
            "germline_distance",
        ],
    threads: 1
    log:
        notebook=rules.windows_report.log.notebook.replace("report", "tune"),
    notebook:
        "../scripts/tune.py.ipynb"


rule windows:
    input:
        expand(
            rules.model_windows.output,
            reads="filtered",
            genome=config["genome"]["name"],
            outdir=config["outdir"],
        ),
