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


rule vcf2bed:
    input:
        vcf=lambda wc: config["genome"][wc.vcf],
        meta=config["donors"],
    output:
        "{outdir}/results/{genome}/vcf2bed/{donor}/{vcf}.bed",
    log:
        "{outdir}/results/{genome}/vcf2bed/{donor}/{vcf}.log",
    conda:
        "../envs/ref.yml"
    shell:
        """
        exec &> {log}

        LIBDID=$(csvcut -t -c "donor_id","libd_id" {input.meta} | csvgrep -c 1 -r "^{wildcards.donor}$" | csvcut -c 2 | tail -n +2)
        # check if the vcf is a xtea vcf
        if [ "{wildcards.vcf}" == "xtea" ]; then
            bcftools view -s "$LIBDID.md" {input.vcf} | \
                bcftools view -i "AC>0" | \
                bcftools query -f '%CHROM\t%POS\t%END\t%FILTER\t%SUBTYPE\t%STRAND\n' -e 'INFO/SUBTYPE ~ "transduction"' | \
                awk -v OFS='\t' '{{print $1, $2-1, $3, $4, $5, $6}}' > {output}

        else
            bcftools view -s "$LIBDID" -i 'INFO/SVTYPE="LINE/L1"' {input.vcf} | \
                bcftools view -i "AC>0" | \
                bcftools query -f '%CHROM\t%0START\t%0END\t%FILTER\t%MEI\t%MESTRAND\n' > {output}
        fi
        """


def get_vcfs(wildcards):
    return {
        "xtea": expand(
            rules.vcf2bed.output,
            vcf="xtea",
            allow_missing=True,
        ),
        "megane_percentile": expand(
            rules.vcf2bed.output,
            vcf="megane_percentile",
            allow_missing=True,
        ),
        "megane_gaussian": expand(
            rules.vcf2bed.output,
            vcf="megane_gaussian",
            allow_missing=True,
        ),
    }


rule coverage:
    input:
        unpack(get_vcfs),
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        fai=config["genome"]["fai"],
        primer_sites=rules.blast_primers.output.bed,
        rmsk=rules.run_rmsk.output.bed,
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

        # filter duplicates and read2
        samtools view -F 1024 -F 128 -F 256 -hb {input.bam} > $tmp_bam
        samtools index $tmp_bam

        echo "Running bedtools coverage for primer sites"
        bedtools coverage -a {input.primer_sites} -b $tmp_bam > {output.primer_sites}

        echo "Running bedtools coverage for xtea"
        bedtools slop -l 100 -r 100 -s -g {input.fai} -i {input.xtea} > $tmp_bed
        bedtools coverage -a $tmp_bed -b $tmp_bam > {output.xtea}

        echo "Running bedtools coverage for megane gaussian"
        bedtools slop -l 100 -r 100 -s -g {input.fai} -i {input.megane_gaussian} > $tmp_bed
        bedtools coverage -a $tmp_bed -b $tmp_bam > {output.megane_gaussian}

        echo "Running bedtools coverage for megane percentile"
        bedtools slop -l 100 -r 100 -s -g {input.fai} -i {input.megane_percentile} > $tmp_bed
        bedtools coverage -a $tmp_bed -b $tmp_bam > {output.megane_percentile}

        echo "Running bedtools coverage for rmsk"
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.rmsk} > $tmp_bed
        bedtools coverage -a $tmp_bed -b $tmp_bam > {output.rmsk}
        """


rule macs3:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
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
        macs3 predictd -i {input.bam} -g hs -f BAMPE 2> {log}
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


# TODO: add closest and overlap labels
rule label:
    input:
        unpack(get_vcfs),
        bulk_regions=lambda wc: expand(
            rules.make_regions.output.bed,
            sample=get_donor_bulk(wc.donor),
            allow_missing=True,
        ),
        cell_regions=rules.make_regions.output.bed,
        primer_sites=rules.blast_primers.output.bed,
        fai=config["genome"]["fai"],
        rmsk=rules.run_rmsk.output.bed,
    output:
        bulk="{outdir}/results/{genome}/{params}/{donor}/{sample}.bulk.bed",
        xtea="{outdir}/results/{genome}/{params}/{donor}/{sample}.xtea.bed",
        megane_gaussian="{outdir}/results/{genome}/{params}/{donor}/{sample}.megane_gaussian.bed",
        megane_percentile="{outdir}/results/{genome}/{params}/{donor}/{sample}.megane_percentile.bed",
        rmsk="{outdir}/results/{genome}/{params}/{donor}/{sample}.rmsk.bed",
        primer_sites="{outdir}/results/{genome}/{params}/{donor}/{sample}.primer_sites.bed",
    conda:
        "../envs/ref.yml"
    log:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}_label.log",
    shell:
        """
        exec &> {log}

        bedtools intersect -a {input.cell_regions} -b {input.bulk_regions} -wa > {output.bulk}

        bedtools intersect -a {input.cell_regions} -b {input.primer_sites} -wa > {output.primer_sites}

        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.xtea} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.xtea}

        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_gaussian} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.megane_gaussian}

        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_percentile} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.megane_percentile}

        # TODO: move this code to rmsk rule filter rmsk for repEnd > 860
        awk '$13 > 860' {input.rmsk} | bedtools slop -l -500 -r 100 -s -g {input.fai} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.rmsk}
        """


rule final_features:
    input:
        regions=rules.make_regions.output.pqt,
        annotations=expand(
            "{outdir}/results/{genome}/{params}/{donor}/{sample}.{annotation}.bed",
            annotation=[
                "bulk",
                "xtea",
                "megane_gaussian",
                "megane_percentile",
                "rmsk",
                "primer_sites",
            ],
            allow_missing=True,
        ),
        en_motif_pos=rules.en_motif.output.pos,
        en_motif_neg=rules.en_motif.output.neg,
    output:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}_labelled.pqt",
    log:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}_labelled.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/final_features.py"


rule eval_regions:
    input:
        cell_regions=rules.final_features.output,
        bulk_regions=lambda wc: expand(
            rules.final_features.output,
            sample=get_donor_bulk(wc.donor),
            allow_missing=True,
        ),
        cell_coverage=expand(
            "{outdir}/results/{genome}/coverage/{donor}/{sample}.{annotation}.bed",
            annotation=[
                "xtea",
                "megane_gaussian",
                "megane_percentile",
                "rmsk",
                "primer_sites",
            ],
            allow_missing=True,
        ),
        bulk_coverage=expand(
            "{outdir}/results/{genome}/coverage/{donor}/{sample}.{annotation}.bed",
            annotation=[
                "xtea",
                "megane_gaussian",
                "megane_percentile",
                "rmsk",
                "primer_sites",
            ],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}_eval.tsv",
    log:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}_eval.log",
    conda:
        "../envs/model.yml"
    script:
        "../scripts/eval_regions.py"


def aggregate(wildcards):
    out = defaultdict(list)
    for d in donors["donor_id"].unique():
        r = expand(
            rules.final_features.output,
            sample=get_donor_cells(d),
            donor=d,
            allow_missing=True,
        )
        out["regions"].extend(r)
        r = expand(
            rules.eval_regions.output,
            sample=get_donor_cells(d),
            donor=d,
            allow_missing=True,
        )
        out["eval"].extend(r)
    return out


rule tune:
    input:
        unpack(aggregate),
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
