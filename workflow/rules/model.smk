rule profile_regions:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
    output:
        pqt="{outdir}/results/{genome}/profile_regions/{donor}/{sample}.pqt",
        stats="{outdir}/results/{genome}/profile_regions/{donor}/{sample}_stats.tsv",
        bed="{outdir}/results/{genome}/profile_regions/{donor}/{sample}.bed",
    log:
        "{outdir}/results/{genome}/profile_regions/{donor}/{sample}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/profile_regions.py"


rule make_regions:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
    output:
        pqt="{outdir}/results/{genome}/{params}/{donor}/{sample}.pqt",
        bed="{outdir}/results/{genome}/{params}/{donor}/{sample}.bed",
    log:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/make_regions.py"


rule coverage:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        fai=config["genome"]["fai"],
        primer_sites=rules.blast_primers.output.bed,
        xtea_vcf=config["genome"]["xtea"],
        megane_percentile_vcf=config["genome"]["megane_percentile"],
        megane_gaussian_vcf=config["genome"]["megane_gaussian"],
        rmsk=rules.run_rmsk.output.bed,
        meta=config["donors"],
    output:
        xtea="{outdir}/results/{genome}/coverage/{donor}/{sample}.xtea.bed",
        megane_gaussian="{outdir}/results/{genome}/coverage/{donor}/{sample}.megane_gaussian.bed",
        megane_percentile="{outdir}/results/{genome}/coverage/{donor}/{sample}.megane_percentile.bed",
        rmsk="{outdir}/results/{genome}/coverage/{donor}/{sample}.rmsk.bed",
        primer_sites="{outdir}/results/{genome}/coverage/{donor}/{sample}.primer_sites.bed",
    log:
        "{outdir}/results/{genome}/coverage/{donor}/{sample}.log",
    conda:
        "../envs/ref.yml"
    shell:
        """
        # send all stdout and stderr to log file
        exec &> {log}

        # setup temp file
        tmp_bed=$(mktemp --suffix .bed)
        tmp_bam=$(mktemp --suffix .bam)
        tmp_bai=$tmp_bam.bai
        trap "rm -f $tmp_bed $tmp_bam $tmp_bai" EXIT

        # filter duplicates
        samtools view -F 1024 -hb {input.bam} > $tmp_bam
        samtools index $tmp_bam

        # grep for wildcards.donor in column titled "donor_id" and get the value in column titled "libd_id"
        LIBDID=$(csvcut -t -c "donor_id","libd_id" {input.meta} | csvgrep -c 1 -m "{wildcards.donor}" | csvcut -c 2 | tail -n +2)
        echo "LIBDID: $LIBDID"

        echo "Running bedtools coverage for primer sites"
        bedtools coverage -a {input.primer_sites} -b $tmp_bam > {output.primer_sites}

        echo "Running bedtools coverage for xtea"
        bcftools view -s $LIBDID.md {input.xtea_vcf} | \
            bcftools view -i "AC>0" | \
            vcf2bed > $tmp_bed
        bedtools coverage -a $tmp_bed -b $tmp_bam > {output.xtea}

        echo "Running bedtools coverage for megane gaussian"
        bcftools view -s $LIBDID -i 'INFO/SVTYPE="LINE/L1"' {input.megane_gaussian_vcf} | \
            bcftools view -i "AC>0" | \
            vcf2bed | \
            bedtools slop -b 200 -g {input.fai} > $tmp_bed
        bedtools coverage -a $tmp_bed -b $tmp_bam > {output.megane_gaussian}

        echo "Running bedtools coverage for megane percentile"
        bcftools view -s $LIBDID -i 'INFO/SVTYPE="LINE/L1"' {input.megane_percentile_vcf} | \
            bcftools view -i "AC>0" | \
            vcf2bed | \
            bedtools slop -b 200 -g {input.fai} > $tmp_bed
        bedtools coverage -a $tmp_bed -b $tmp_bam > {output.megane_percentile}

        echo "Running bedtools coverage for rmsk"
        bedtools slop -l 0 -r 200 -s -g {input.fai} -i {input.rmsk} > $tmp_bed
        bedtools coverage -a $tmp_bed -b $tmp_bam > {output.rmsk}
        """


rule macs3:
    input:
        bam=rules.rmdup.output.bam,
        bai=rules.rmdup.output.idx,
        dup_bam=rules.sort.output[0],
        dup_bai=rules.index.output[0],
    output:
        xls="{outdir}/results/{genome}/macs3/{donor}/{sample}_peaks.xls",
        narrowPeak="{outdir}/results/{genome}/macs3/{donor}/{sample}_peaks.narrowPeak",
        summits="{outdir}/results/{genome}/macs3/{donor}/{sample}_summits.bed",
        cutoff_analysis="{outdir}/results/{genome}/macs3/{donor}/{sample}_cutoff_analysis.txt",
    log:
        "{outdir}/results/{genome}/macs3/{donor}/{sample}.log",
    conda:
        "../envs/features.yml"
    shell:
        """
        macs3 predictd -i {input.dup_bam} -g hs -f BAMPE 2> {log}
        insert_size=$(grep "Average insertion length" {log} | cut -d " " -f 18)
        macs3 callpeak \
            -t {input.bam} -f BAM -g hs \
            --nomodel \
            --extsize $insert_size \
            --shift 0 \
            --llocal 1000 --slocal 1000 \
            -p 0.05 \
            --outdir $(dirname {output[0]}) \
            -n {wildcards.sample} \
            --cutoff-analysis \
            --keep-dup all 2>> {log}
        """


with open("resources/bad_cells.txt", "r") as f:
    bad_cells = [line.strip() for line in f.readlines()]


def get_donor_cells(donor):
    cells = samples.loc[
        (samples["donor_id"] == donor) & (~samples["sample_id"].str.contains("gDNA"))
    ]["sample_id"].values
    cells = [c for c in cells if c not in bad_cells]
    return cells


def get_donor_bulk(donor):
    return samples.loc[
        (samples["donor_id"] == donor) & (samples["sample_id"].str.contains("gDNA"))
    ]["sample_id"].values


def get_label_input(wildcards):
    return {
        "bulk": expand(
            rules.make_regions.output.pqt,
            sample=get_donor_bulk(wildcards.donor),
            allow_missing=True,
        ),
        "cells": expand(
            rules.make_regions.output.pqt,
            sample=get_donor_cells(wildcards.donor),
            allow_missing=True,
        ),
    }


rule label:
    input:
        unpack(get_label_input),
        en_motif_pos=rules.en_motif.output.pos,
        en_motif_neg=rules.en_motif.output.neg,
        primer_sites=rules.blast_primers.output.bed,
        xtea_vcf=config["genome"]["xtea"],
        megane_percentile_vcf=config["genome"]["megane_percentile"],
        megane_gaussian_vcf=config["genome"]["megane_gaussian"],
        rmsk=rules.run_rmsk.output.bed,
    output:
        cells="{outdir}/results/{genome}/{params}/{donor}.pqt",
        bulk="{outdir}/results/{genome}/{params}/{donor}_bulk.pqt",
    log:
        "{outdir}/results/{genome}/{params}/{donor}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/label.py"


def get_eval_regions_input(wildcards):
    out = {}
    out["cells_regions"] = rules.label.output.cells
    for cov in ["xtea", "megane_gaussian", "megane_percentile", "rmsk", "primer_sites"]:
        out[f"cells_coverage_{cov}"] = expand(
            rules.coverage.output[cov],
            sample=get_donor_cells(wildcards.donor),
            allow_missing=True,
        )
        out[f"bulk_coverage_{cov}"] = expand(
            rules.coverage.output[cov],
            sample=get_donor_bulk(wildcards.donor),
            allow_missing=True,
        )
    return {
        "bulk_regions": rules.label.output.bulk,
        "cells_regions": rules.label.output.cells,
        "cells_coverage": expand(
            rules.coverage.output,
            sample=get_donor_cells(wildcards.donor),
            allow_missing=True,
        ),
        "bulk_coverage": expand(
            rules.coverage.output,
            sample=get_donor_bulk(wildcards.donor),
            allow_missing=True,
        ),
    }


rule eval_regions:
    input:
        unpack(get_eval_regions_input),
    output:
        "{outdir}/results/{genome}/{params}/{donor}_eval.pdf",
    log:
        "{outdir}/results/{genome}/{params}/{donor}_eval.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/eval_regions.py"


rule tune:
    input:
        expand(
            rules.label.output,
            donor=donors["donor_id"].unique(),
            allow_missing=True,
        ),
    output:
        model="{outdir}/results/{genome}/{params}/model/model.pkl",
        best_hp="{outdir}/results/{genome}/{params}/model/best_hp.json",
        history="{outdir}/results/{genome}/{params}/model/history.log",
        tuning_curve="{outdir}/results/{genome}/{params}/model/tuning_curve.png",
        precision_recall_curve="{outdir}/results/{genome}/{params}/model/precision_recall_curve.png",
    log:
        "{outdir}/results/{genome}/{params}/tune.log",
    conda:
        "../envs/model.yml"
    params:
        max_iter=50,
    threads: 32
    script:
        "../scripts/tune.py"
