rule wgs_bp_reads:
    input:
        bam="resources/{genome}/single_cell_wgs/1/{s}/{s}.md.bam",
        breakpoints="resources/{genome}/single_cell_wgs/1/{s}/breakpoint_pairs.txt.gz",
    output:
        reads=temp("{outdir}/results/{genome}/single_cell_wgs/1/{s}/{s}.bp_reads.txt"),
        bam="{outdir}/results/{genome}/single_cell_wgs/1/{s}/{s}.bp_reads.bam",
        bai="{outdir}/results/{genome}/single_cell_wgs/1/{s}/{s}.bp_reads.bam.bai",
    log:
        "{outdir}/results/{genome}/single_cell_wgs/1/{s}/{s}.log",
    conda:
        "../envs/align.lock.yml"
    threads: 8
    shell:
        """
        touch {log} & exec 1> {log} 2>&1
        zcat {input.breakpoints} | \
            grep "LINE/L1" | \
            awk '{{print $9;$10}}' | \
            tr ';' '\n' | \
            awk -F '/' '{{print $3}}' > {output.reads}
        samtools view -hb -N {output.reads} -@ {threads} {input.bam} > {output.bam}
        samtools index {output.bam}
        """


rule rsync_to_brainome:
    input:
        single_cell_wgs_bams=[
            "resources/{genome}/single_cell_wgs/1/A8/A8.md.bam",
            "resources/{genome}/single_cell_wgs/1/B3/B3.md.bam",
            "resources/{genome}/single_cell_wgs/1/D6/D6.md.bam",
        ],
        slav_bams=expand(
            rules.sort.output[0],
            donor="1",
            sample=["ush1_A8_S178", "ush1_B3_S140", "ush1_D6_S165", "gDNA_usd1"],
            allow_missing=True,
        ),
        single_cell_wgs_bp_bams=expand(
            rules.wgs_bp_reads.output.bam, s=["A8", "B3", "D6"], allow_missing=True
        ),
        bulk_bam="resources/{genome}/wgs_calls/30x/LIBD73/LIBD73.md.bam",
    output:
        "{outdir}/results/{genome}/{reads}/peaks/rsync_to_brainome.done",
    shell:
        """
        for f in {input.single_cell_wgs_bams} {input.single_cell_wgs_bp_bams}; do
            rsync -avzP $f brainome:/mysqlpool/mcuoco/for_igv/single_cell_30x_wgs/$(basename $f)
            rsync -avzP $f.bai brainome:/mysqlpool/mcuoco/for_igv/single_cell_30x_wgs/$(basename $f).bai
        done

        for f in {input.slav_bams}; do
            rsync -avzP $f brainome:/mysqlpool/mcuoco/for_igv/slavseq/$(basename $f)
            rsync -avzP $f.bai brainome:/mysqlpool/mcuoco/for_igv/slavseq/$(basename $f).bai
        done

        rsync -avzP {input.bulk_bam} brainome:/mysqlpool/mcuoco/for_igv/bulk_30x_wgs_calls/1/$(basename {input.bulk_bam})
        rsync -avzP {input.bulk_bam}.bai brainome:/mysqlpool/mcuoco/for_igv/bulk_30x_wgs_calls/1/$(basename {input.bulk_bam}).bai

        touch {output}
        """


rule wgs_igv:
    input:
        rsync=rules.rsync_to_brainome.output,
        my_predictions=expand(
            rules.tune.output.predictions,
            features=feature_sets.keys(),
            allow_missing=True,
        ),
        apua_predictions="resources/meiyan_apua_validation/U01_L1_primers.csv",
        shap_values=expand(
            rules.tune.output.shap_values,
            features=feature_sets.keys(),
            allow_missing=True,
        ),
        wgs_bams=[
            "resources/{genome}/single_cell_wgs/1/A8/A8.md.bam",
            "resources/{genome}/single_cell_wgs/1/B3/B3.md.bam",
            "resources/{genome}/single_cell_wgs/1/D6/D6.md.bam",
            "resources/{genome}/wgs_calls/30x/LIBD73/LIBD73.md.bam",
        ],
        slav_bams=expand(
            rules.sort.output[0],
            donor="1",
            sample=["ush1_A8_S178", "ush1_B3_S140", "ush1_D6_S165", "gDNA_usd1"],
            allow_missing=True,
        ),
        breakpoints=[
            "resources/{genome}/single_cell_wgs/1/A8/breakpoint_pairs_pooled_all.txt.gz",
            "resources/{genome}/single_cell_wgs/1/B3/breakpoint_pairs_pooled_all.txt.gz",
            "resources/{genome}/single_cell_wgs/1/D6/breakpoint_pairs_pooled_all.txt.gz",
            "resources/{genome}/wgs_calls/30x/LIBD74/breakpoint_pairs_pooled_all.txt.gz",
        ],
        wgs_bp_bams=expand(
            rules.wgs_bp_reads.output.bam, s=["A8", "B3", "D6"], allow_missing=True
        ),
    output:
        "{outdir}/results/{genome}/{reads}/peaks/igv.bat",
    log:
        notebook="{outdir}/results/{genome}/{reads}/peaks/wgs_igv.ipynb",
    notebook:
        "../scripts/wgs_igv.py.ipynb"
