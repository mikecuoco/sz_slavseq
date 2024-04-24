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
        pqt="{outdir}/results/{genome}/peaks/{donor}/{sample}.pqt",
        bed="{outdir}/results/{genome}/peaks/{donor}/{sample}.bed",
    log:
        "{outdir}/results/{genome}/peaks/{donor}/{sample}.log",
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


# TODO: add closest and overlap labels
rule label:
    input:
        unpack(get_vcfs),
        bulk_regions=lambda wc: expand(
            rules.make_regions.output.bed,
            sample=bulk.loc[wc.donor, "sample_id"],
            allow_missing=True,
        ),
        cell_regions=rules.make_regions.output.bed,
        primer_sites=rules.blast_primers.output.bed,
        fai=config["genome"]["fai"],
        l1hs=rules.filter_rmsk.output.l1hs,
        l1pa2=rules.filter_rmsk.output.l1pa2,
        l1pa3=rules.filter_rmsk.output.l1pa3,
        l1pa4=rules.filter_rmsk.output.l1pa4,
        l1pa5=rules.filter_rmsk.output.l1pa5,
        l1pa6=rules.filter_rmsk.output.l1pa6,
    output:
        bulk="{outdir}/results/{genome}/peaks/{donor}/{sample}.bulk.bed",
        xtea="{outdir}/results/{genome}/peaks/{donor}/{sample}.xtea.bed",
        megane_gaussian="{outdir}/results/{genome}/peaks/{donor}/{sample}.megane_gaussian.bed",
        megane_percentile="{outdir}/results/{genome}/peaks/{donor}/{sample}.megane_percentile.bed",
        primer_sites="{outdir}/results/{genome}/peaks/{donor}/{sample}.primer_sites.bed",
        l1hs="{outdir}/results/{genome}/peaks/{donor}/{sample}.l1hs.bed",
        l1pa2="{outdir}/results/{genome}/peaks/{donor}/{sample}.l1pa2.bed",
        l1pa3="{outdir}/results/{genome}/peaks/{donor}/{sample}.l1pa3.bed",
        l1pa4="{outdir}/results/{genome}/peaks/{donor}/{sample}.l1pa4.bed",
        l1pa5="{outdir}/results/{genome}/peaks/{donor}/{sample}.l1pa5.bed",
        l1pa6="{outdir}/results/{genome}/peaks/{donor}/{sample}.l1pa6.bed",
    conda:
        "../envs/ref.yml"
    log:
        "{outdir}/results/{genome}/peaks/{donor}/{sample}_label.log",
    shell:
        """
        exec &> {log}

        bedtools intersect -a {input.cell_regions} -b {input.bulk_regions} -wa > {output.bulk}
        bedtools intersect -a {input.cell_regions} -b {input.primer_sites} -wa > {output.primer_sites}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.xtea} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.xtea}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_gaussian} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.megane_gaussian}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_percentile} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.megane_percentile}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1hs} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.l1hs}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa2} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.l1pa2}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa3} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.l1pa3}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa4} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.l1pa4}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa5} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.l1pa5}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa6} | bedtools intersect -a {input.cell_regions} -b - -wa > {output.l1pa6}
        """


rule final_features:
    input:
        regions=rules.make_regions.output.pqt,
        annotations=expand(
            "{outdir}/results/{genome}/peaks/{donor}/{sample}.{annotation}.bed",
            annotation=[
                "bulk",
                "xtea",
                "megane_gaussian",
                "megane_percentile",
                "primer_sites",
                "l1hs",
                "l1pa2",
                "l1pa3",
                "l1pa4",
                "l1pa5",
                "l1pa6",
            ],
            allow_missing=True,
        ),
        en_motif_pos=rules.en_motif.output.pos,
        en_motif_neg=rules.en_motif.output.neg,
    output:
        "{outdir}/results/{genome}/peaks/{donor}/{sample}_labelled.pqt",
    log:
        "{outdir}/results/{genome}/peaks/{donor}/{sample}_labelled.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/final_features.py"


rule eval_regions:
    input:
        cell_regions=rules.final_features.output,
        bulk_regions=lambda wc: expand(
            rules.final_features.output,
            sample=bulk.loc[wc.donor, "sample_id"],
            allow_missing=True,
        ),
        cell_coverage=expand(
            "{outdir}/results/{genome}/coverage/{donor}/{sample}.{annotation}.bed",
            annotation=[
                "xtea",
                "megane_gaussian",
                "megane_percentile",
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
                "primer_sites",
            ],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/{genome}/peaks/{donor}/{sample}_eval.tsv",
    log:
        "{outdir}/results/{genome}/peaks/{donor}/{sample}_eval.log",
    conda:
        "../envs/model.yml"
    script:
        "../scripts/eval_regions.py"


rule peak_coverage:
    input:
        cell=rules.final_features.output,
        bulk=lambda wc: expand(
            rules.final_features.output,
            sample=bulk.loc[wc.donor, "sample_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/{genome}/peaks/{donor}/{sample}_coverage.pqt",
    log:
        "{outdir}/results/{genome}/peaks/{donor}/{sample}_coverage.log",
    conda:
        "../envs/model.yml"
    script:
        "../scripts/peak_coverage.py"


rule germline_report:
    input:
        bulk=expand(
            rules.final_features.output,
            zip,
            sample=bulk["sample_id"],
            donor=bulk["donor_id"],
            allow_missing=True,
        ),
        peak_coverage=expand(
            rules.peak_coverage.output,
            zip,
            donor=bulk["donor_id"],
            sample=bulk["sample_id"],
            allow_missing=True,
        ),
        megane=expand(
            rules.vcf2bed.output,
            vcf="megane_gaussian",
            donor=donors["donor_id"].unique(),
            allow_missing=True,
        ),
        flagstat=expand(
            rules.flagstat.output,
            zip,
            donor=bulk["donor_id"],
            sample=bulk["sample_id"],
            allow_missing=True,
        ),
        meta=config["donors"],
    output:
        "{outdir}/results/{genome}/peaks/bulk.pqt",
    conda:
        "../envs/model.yml"
    log:
        notebook="{outdir}/results/{genome}/peaks/germline_report.ipynb",
    params:
        pos_label="megane_gaussian",
    notebook:
        "../scripts/germline_report.py.ipynb"


rule regions_report:
    input:
        cells=expand(
            rules.final_features.output,
            zip,
            donor=cells["donor_id"],
            sample=cells["sample_id"],
            allow_missing=True,
        ),
        peak_coverage=expand(
            rules.peak_coverage.output,
            zip,
            donor=cells["donor_id"],
            sample=cells["sample_id"],
            allow_missing=True,
        ),
        bulk=rules.germline_report.output,
        flagstat=expand(
            rules.flagstat.output,
            zip,
            donor=cells["donor_id"],
            sample=cells["sample_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/{genome}/peaks/data.pqt",
    conda:
        "../envs/model.yml"
    log:
        notebook="{outdir}/results/{genome}/peaks/regions_report.ipynb",
    params:
        pos_label="megane_gaussian",
    notebook:
        "../scripts/regions_report.py.ipynb"


rule tune:
    input:
        data=rules.regions_report.output,
    output:
        model="{outdir}/results/{genome}/peaks/model/model.pkl",
        best_hp="{outdir}/results/{genome}/peaks/model/best_hp.json",
        history="{outdir}/results/{genome}/peaks/model/history.log",
    log:
        "{outdir}/results/{genome}/peaks/model/tune.log",
    conda:
        "../envs/model.yml"
    params:
        max_iter=100,
        random_state=1,
    threads: 32
    script:
        "../scripts/tune.py"


rule model_report:
    input:
        data=rules.regions_report.output,
        bulk=rules.germline_report.output,
        history=rules.tune.output.history,
        best_hp=rules.tune.output.best_hp,
        model=rules.tune.output.model,
    output:
        "{outdir}/results/{genome}/peaks/model/predictions.pqt",
    conda:
        "../envs/model.yml"
    params:
        random_state=1,
    log:
        notebook="{outdir}/results/{genome}/peaks/model_report.ipynb",
    notebook:
        "../scripts/model_report.py.ipynb"


rule calls_report:
    input:
        data=rules.model_report.output,
        bulk=rules.germline_report.output,
    output:
        expand(
            "{outdir}/results/{genome}/peaks/{donor}/calls.bed",
            donor=donors["donor_id"].unique(),
            allow_missing=True,
        ),
    conda:
        "../envs/model.yml"
    log:
        notebook="{outdir}/results/{genome}/peaks/calls_report.ipynb",
    notebook:
        "../scripts/calls_report.py.ipynb"


def get_donor_bams(wildcards):
    return expand(
        rules.sort.output,
        sample=samples.loc[wildcards.donor, "sample_id"],
        donor=wildcards.donor,
        allow_missing=True,
    )


rule make_igv_batch_script:
    input:
        regions="{outdir}/results/{genome}/peaks/{donor}/calls.bed",
        bams=get_donor_bams,
        rmsk=rules.run_rmsk.output.bed,
        megane_gaussian="{outdir}/results/{genome}/vcf2bed/{donor}/megane_gaussian.bed",
        primer_sites=rules.blast_primers.output.bed,
        fasta=config["genome"]["fasta"],
        fai=config["genome"]["fai"],
    output:
        "{outdir}/results/{genome}/peaks/{donor}/igv_snapshots.bat",
    log:
        "{outdir}/results/{genome}/peaks/{donor}/make_igv_batch_script.log",
    conda:
        "../envs/igv.yml"
    params:
        maxPanelHeight=200,
        colorBy="READ_STRAND",
    script:
        "../scripts/igv_snapshots.py"


rule igv_download:
    output:
        "resources/IGV_Linux_2.17.4/igv.sh",
    params:
        runtime="600",
        memory="1G",
    shell:
        """
        cd resources
        curl https://data.broadinstitute.org/igv/projects/downloads/2.17/IGV_Linux_2.17.4_WithJava.zip > IGV.zip; unzip IGV.zip
        """


rule igv_snapshots:
    input:
        script=rules.make_igv_batch_script.output,
        igv=rules.igv_download.output,
    output:
        directory("{outdir}/results/{genome}/peaks/{donor}/snapshots"),
    log:
        "{outdir}/results/{genome}/peaks/{donor}/igv_snapshots.log",
    shell:
        """
        xvfb-run --auto-servernum {input.igv} -b {input.script} &> {log}
        """
