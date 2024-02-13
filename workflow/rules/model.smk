rule make_regions:
    input:
        bam=rules.sambamba_sort.output[0],
        bai=rules.sambamba_index.output[0],
    output:
        pqt="{outdir}/results/{genome}/{params}/{donor}/{sample}.pqt",
    log:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/make_regions.py"


rule coverage:
    input:
        bam=rules.rmdup.output.bam,
        bai=rules.rmdup.output.idx,
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
        trap "rm -f $tmp_bed" EXIT

        # grep for wildcards.donor in first column, return 5th column
        LIBDID=$(grep ^{wildcards.donor} {input.meta} | cut -f 5)
        echo "LIBDID: $LIBDID"

        echo "Running bedtools coverage for primer sites"
        bedtools coverage -a {input.primer_sites} -b {input.bam} > {output.primer_sites}

        echo "Running bedtools coverage for xtea"
        bcftools view -s $LIBDID.md {input.xtea_vcf} | \
            bcftools view -i "AC>0" | \
            vcf2bed > $tmp_bed
        bedtools coverage -a $tmp_bed -b {input.bam} > {output.xtea}

        echo "Running bedtools coverage for megane gaussian"
        bcftools view -s $LIBDID -i 'INFO/SVTYPE="LINE/L1"' {input.megane_gaussian_vcf} | \
            bcftools view -i "AC>0" | \
            vcf2bed | \
            bedtools slop -b 200 -g {input.fai} > $tmp_bed
        bedtools coverage -a $tmp_bed -b {input.bam} > {output.megane_gaussian}

        echo "Running bedtools coverage for megane percentile"
        bcftools view -s $LIBDID -i 'INFO/SVTYPE="LINE/L1"' {input.megane_percentile_vcf} | \
            bcftools view -i "AC>0" | \
            vcf2bed | \
            bedtools slop -b 200 -g {input.fai} > $tmp_bed
        bedtools coverage -a $tmp_bed -b {input.bam} > {output.megane_percentile}

        echo "Running bedtools coverage for rmsk"
        bedtools slop -l 0 -r 200 -s -g {input.fai} -i {input.rmsk} > $tmp_bed
        bedtools coverage -a $tmp_bed -b {input.bam} > {output.rmsk}
        """


rule macs3:
    input:
        bam=rules.rmdup.output.bam,
        bai=rules.rmdup.output.idx,
        dup_bam=rules.sambamba_sort.output[0],
        dup_bai=rules.sambamba_index.output[0],
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


def get_donor_cells(wildcards):
    cells = samples.loc[
        (samples["donor_id"] == wildcards.donor)
        & (~samples["sample_id"].str.contains("gDNA"))
    ]["sample_id"].values
    cells = [c for c in cells if c not in bad_cells]

    return expand(
        rules.make_regions.output.pqt,
        # rules.macs3.output.narrowPeak,
        sample=cells,
        allow_missing=True,
    )


rule profile_regions:
    input:
        bam=rules.sambamba_sort.output[0],
        bai=rules.sambamba_index.output[0],
    output:
        pqt="{outdir}/results/{genome}/profile_regions/{donor}/{sample}.pqt",
        stats="{outdir}/results/{genome}/profile_regions/{donor}/{sample}_stats.tsv",
    log:
        "{outdir}/results/{genome}/profile_regions/{donor}/{sample}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/profile_regions.py"


def get_donor_bulk(wildcards):
    bulk = samples.loc[
        (samples["donor_id"] == wildcards.donor)
        & (samples["sample_id"].str.contains("gDNA"))
    ]["sample_id"].values
    bulk = [c for c in bulk if c not in bad_cells]
    return expand(
        rules.make_regions.output.pqt,
        # rules.coverage.output.xtea,
        # rules.macs3.output.narrowPeak,
        sample=bulk,
        allow_missing=True,
    )


rule label:
    input:
        cells=get_donor_cells,
        bulk=get_donor_bulk,
        en_motif_pos=rules.en_motif.output.pos,
        en_motif_neg=rules.en_motif.output.neg,
        primer_sites=rules.blast_primers.output.bed,
        xtea_vcf=config["genome"]["xtea"],
        megane_percentile_vcf=config["genome"]["megane_percentile"],
        megane_gaussian_vcf=config["genome"]["megane_gaussian"],
        rmsk=rules.run_rmsk.output.bed,
    output:
        "{outdir}/results/{genome}/{params}/{donor}.pqt",
    log:
        "{outdir}/results/{genome}/{params}/{donor}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/label.py"


# TODO: add rule to generate report on regions
# By cells, donors, and labels
# 1. # regions
# 2. if peaks, length
# 3. reads / region


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
