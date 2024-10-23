rule make_breakpoints:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        chrom_sizes=config["genome"]["genome"],
    output:
        bed="{outdir}/results/{genome}/{reads}/breakpoints/{donor}/{sample}.bed.gz",
        tbi="{outdir}/results/{genome}/{reads}/breakpoints/{donor}/{sample}.bed.gz.tbi",
    log:
        "{outdir}/results/{genome}/{reads}/breakpoints/{donor}/{sample}.make_breakpoints.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/make_breakpoints.py"


def get_breakpoints(wildcards):
    if wildcards.group == "bulk":
        return expand(
            rules.make_breakpoints.output.bed,
            zip,
            sample=bulk["sample_id"],
            donor=bulk["donor_id"],
            allow_missing=True,
        )
    elif wildcards.group == "cells":
        return expand(
            rules.make_breakpoints.output.bed,
            zip,
            sample=cells["sample_id"],
            donor=cells["donor_id"],
            allow_missing=True,
        )
    else:
        raise ValueError(f"Unknown breakpoint group: {wildcards.group}")


rule merge_breakpoints:
    input:
        get_breakpoints,
    output:
        bed="{outdir}/results/{genome}/{reads}/breakpoints/{group}.bed.gz",
        tbi="{outdir}/results/{genome}/{reads}/breakpoints/{group}.bed.gz.tbi",
    log:
        "{outdir}/results/{genome}/{reads}/breakpoints/{group}_merge.log",
    conda:
        "../envs/align.lock.yml"
    shell:
        """
        exec &>> {log}

        # write header
        bed={output.bed}
        bed=${{bed%.gz}}
        echo "bed is $bed"
        (zcat {input[0]} | head -n 1 || true) > $bed
        echo "merging and sorting..."
        # inputs=({input})
        for i in {input}; do
            zcat $i | tail -n +2
        done | sort -k1,1 -k2,2n >> $bed

        echo "compressing"
        bgzip $bed
        echo "index..."
        tabix -f -p bed {output.bed}
        """


label_out = {}
for label in [
    "bulk",
    "megane",
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
    label_out[label] = rules.merge_breakpoints.output.bed.replace(
        ".bed.gz", f".{label}.bed.gz"
    )


# TODO: add closest and overlap labels
rule label_breakpoints:
    input:
        **rmsk_anno,
        **knrgl_anno,
        bulk_regions="{outdir}/results/{genome}/{reads}/breakpoints/bulk.bed.gz",
        cell_regions="{outdir}/results/{genome}/{reads}/breakpoints/{group}.bed.gz",
        primer_sites=rules.blast_primers.output.bed,
        fai=config["genome"]["fai"],
    output:
        **label_out,
    conda:
        "../envs/ref.lock.yml"
    log:
        rules.merge_breakpoints.log[0].replace("merge", "label"),
    shell:
        """
        exec &>> {log}

        echo "labelling bulk regions"
        bedtools intersect -a {input.bulk_regions} -b {input.cell_regions} -wao | bgzip -ci -I {output.bulk}.gzi > {output.bulk}
        echo "labelling primer sites"
        bedtools intersect -a {input.primer_sites} -b {input.cell_regions} -wao | bgzip -ci -I {output.primer_sites}.gzi > {output.primer_sites}
        echo "labelling megane calls"
        bedtools intersect -a {input.megane} -b {input.cell_regions} -wao | bgzip -ci -I {output.megane}.gzi > {output.megane}
=        echo "labelling reference l1hs"
        bedtools intersect -a {input.l1hs} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1hs}.gzi > {output.l1hs}
        echo "labelling reference l1pa2"
        bedtools intersect -a {input.l1pa2} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa2}.gzi > {output.l1pa2}
        echo "labelling reference l1pa3"
        bedtools intersect -a {input.l1pa3} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa3}.gzi > {output.l1pa3}
        echo "labelling reference l1pa4"
        bedtools intersect -a {input.l1pa4} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa4}.gzi > {output.l1pa4}
        echo "labelling reference l1pa5"
        bedtools intersect -a {input.l1pa5} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa5}.gzi > {output.l1pa5}
        echo "labelling reference l1pa6"
        bedtools intersect -a {input.l1pa6} -b {input.cell_regions} -wao | bgzip -ci -I {output.l1pa6}.gzi > {output.l1pa6}
        echo "labelling reference polyA"
        bedtools intersect -a {input.polyA} -b {input.cell_regions} -wao | bgzip -ci -I {output.polyA}.gzi > {output.polyA}
        echo "labelling reference polyT"
        bedtools intersect -a {input.polyT} -b {input.cell_regions} -wao | bgzip -ci -I {output.polyT}.gzi > {output.polyT}
        """


rule merge_breakpoint_labels:
    input:
        regions="{outdir}/results/{genome}/{reads}/breakpoints/{group}.bed.gz",
        annotations=expand(
            rules.label_breakpoints.output.bulk.replace("bulk", "{annotation}"),
            annotation=label_out,
            allow_missing=True,
        ),
    output:
        rules.label_breakpoints.output.bulk.replace("bulk.bed.gz", "labelled.bed.gz"),
    log:
        rules.label_breakpoints.log[0].replace("label", "merge"),
    conda:
        "../envs/features.yml"
    script:
        "../scripts/merge_labels.py"


rule breakpoints_sensitivity:
    input:
        bulk="{outdir}/results/{genome}/{reads}/breakpoints/bulk.labelled.bed.gz",
        cells="{outdir}/results/{genome}/{reads}/breakpoints/cells.labelled.bed.gz",
        l1hs_rmsk="resources/{genome}/{genome}.fasta.rmsk.l1hs.bed",
        megane=expand(
            rules.vcf2bed.output,
            vcf="megane_gaussian",
            donor=bulk["donor_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/{genome}/{reads}/breakpoints/sensitivity.ipynb",
    log:
        notebook="{outdir}/results/{genome}/{reads}/breakpoints/sensitivity.ipynb",
    conda:
        "../envs/model.yml"
    notebook:
        "../scripts/sensitivity.py.ipynb"


rule bulk_breakpoints_report:
    input:
        bulk="{outdir}/results/{genome}/{reads}/breakpoints/bulk.labelled.bed.gz",
    output:
        "{outdir}/results/{genome}/{reads}/breakpoints/bulk_report.ipynb",
    log:
        notebook="{outdir}/results/{genome}/{reads}/breakpoints/bulk_report.ipynb",
    conda:
        "../envs/model.yml"
    notebook:
        "../scripts/bulk_peaks_report.py.ipynb"


rule breakpoints:
    input:
        expand(
            rules.breakpoints_sensitivity.output,
            reads="filtered",
            genome=config["genome"]["name"],
            outdir=config["outdir"],
        ),
