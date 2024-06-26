rule make_regions:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        en_motif_pos=rules.en_motif.output.pos,
        en_motif_neg=rules.en_motif.output.neg,
    output:
        pqt=rules.bwa_mem.output[0]
        .replace(".genome.bam", ".pqt")
        .replace("align", "peaks"),
        bed=rules.bwa_mem.output[0]
        .replace(".genome.bam", ".bed")
        .replace("align", "peaks"),
    log:
        "{outdir}/results/{genome}/{reads}/peaks/{donor}/{sample}.make_regions.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/make_regions.py"


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
        "megane_breakpoints": expand(
            rules.vcf2bed.output,
            vcf="megane_breakpoints",
            allow_missing=True,
        ),
        "graffite": expand(
            rules.vcf2bed.output,
            vcf="graffite",
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
        "../envs/ref.lock.yml"
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
]:
    label_out[label] = rules.make_regions.output.bed.replace(".bed", f".{label}.bed")


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
        **label_out,
    conda:
        "../envs/ref.lock.yml"
    log:
        rules.make_regions.log[0].replace("make_regions", "label"),
    shell:
        """
        exec &>> {log}

        bedtools intersect -a {input.bulk_regions} -b {input.cell_regions} -wao > {output.bulk}
        bedtools intersect -a {input.primer_sites} -b {input.cell_regions} -wao > {output.primer_sites}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.xtea} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.xtea}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_gaussian} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.megane_gaussian}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_percentile} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.megane_percentile}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.megane_breakpoints} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.megane_breakpoints}
        bedtools slop -l 100 -r 100 -g {input.fai} -i {input.graffite} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.graffite}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1hs} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.l1hs}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa2} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.l1pa2}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa3} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.l1pa3}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa4} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.l1pa4}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa5} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.l1pa5}
        bedtools slop -l -500 -r 100 -s -g {input.fai} -i {input.l1pa6} | bedtools intersect -a - -b {input.cell_regions} -wao > {output.l1pa6}
        """


rule merge_labels:
    input:
        regions=rules.make_regions.output.pqt,
        annotations=expand(
            rules.label.output.bulk.replace("bulk", "{annotation}"),
            annotation=[
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
            ],
            allow_missing=True,
        ),
    output:
        data=rules.label.output.bulk.replace("bulk.bed", "labelled.pqt"),
        coverage=rules.label.output.bulk.replace("bulk.bed", "coverage.tsv"),
    log:
        rules.label.log[0].replace("label", "merge"),
    conda:
        "../envs/features.yml"
    script:
        "../scripts/merge_labels.py"


rule germline_report:
    input:
        bulk=expand(
            rules.merge_labels.output.data,
            zip,
            sample=bulk["sample_id"],
            donor=bulk["donor_id"],
            allow_missing=True,
        ),
        bulk_coverage=expand(
            rules.merge_labels.output.coverage,
            zip,
            donor=bulk["donor_id"],
            sample=bulk["sample_id"],
            allow_missing=True,
        ),
        megane=expand(
            rules.vcf2bed.output,
            vcf=[
                "xtea",
                "megane_gaussian",
                "megane_percentile",
                "megane_breakpoints",
                "graffite",
            ],
            donor=donors["donor_id"].unique(),
            allow_missing=True,
        ),
        flagstat=expand(
            rules.flagstat.output,
            zip,
            sample=bulk["sample_id"],
            donor=bulk["donor_id"],
            allow_missing=True,
        ),
        meta=config["donors"],
    output:
        rules.reads_report.output[0].replace("reads_report.ipynb", "peaks/bulk.pqt"),
    conda:
        "../envs/model.yml"
    log:
        notebook=rules.reads_report.output[0].replace(
            "reads_report.ipynb", "peaks/germline_report.ipynb"
        ),
    notebook:
        "../scripts/germline_report.py.ipynb"


rule multiinter:
    input:
        bed=expand(
            rules.make_regions.output.bed,
            zip,
            sample=samples["sample_id"],
            donor=samples["donor_id"],
            allow_missing=True,
        ),
    output:
        rules.germline_report.output[0].replace("bulk.pqt", "multiinter.bed"),
    log:
        rules.germline_report.output[0].replace("bulk.pqt", "multiinter.log"),
    conda:
        "../envs/ref.lock.yml"
    shell:
        """
        names=()
        for bed in {input.bed}; do
            names+=($(basename $bed .bed))
        done
        echo ${{names[@]}} > {log}
        bedtools multiinter -i {input.bed} -names ${{names[@]}} > {output} 2>> {log}
        """


rule regions_report:
    input:
        cell_peaks=expand(
            rules.merge_labels.output.data,
            zip,
            donor=cells["donor_id"],
            sample=cells["sample_id"],
            allow_missing=True,
        ),
        cell_coverage=expand(
            rules.merge_labels.output.coverage,
            zip,
            sample=cells["sample_id"],
            donor=cells["donor_id"],
            allow_missing=True,
        ),
        bulk_peaks=rules.germline_report.output,
        multiinter=rules.multiinter.output,
    output:
        rules.germline_report.output[0].replace("bulk.pqt", "data.pqt"),
    conda:
        "../envs/model.yml"
    log:
        notebook=rules.germline_report.log.notebook.replace(
            "germline_report", "regions_report"
        ),
    notebook:
        "../scripts/regions_report.py.ipynb"


rule tune:
    input:
        data=rules.regions_report.output,
        cell_coverage=expand(
            rules.merge_labels.output.coverage,
            zip,
            sample=cells["sample_id"],
            donor=cells["donor_id"],
            allow_missing=True,
        ),
    output:
        model=rules.regions_report.output[0].replace("data.pqt", "{features}/model.pkl"),
        best_hp=rules.regions_report.output[0].replace(
            "data.pqt", "{features}/best_hp.json"
        ),
        history=rules.regions_report.output[0].replace(
            "data.pqt", "{features}/history.log"
        ),
        predictions=rules.regions_report.output[0].replace(
            "data.pqt", "{features}/predictions.pqt"
        ),
    conda:
        "../envs/model_gpu.yml" if shutil.which("nvidia-smi") else "../envs/model.yml"
    params:
        max_iter=100,
        random_state=1,
        features=lambda wc: feature_sets[wc.features],
    threads: 1
    log:
        notebook=rules.regions_report.log.notebook.replace(
            "regions_report", "{features}/tune"
        ),
    notebook:
        "../scripts/tune.py.ipynb"


rule model_report:
    input:
        data=rules.regions_report.output,
        bulk=rules.germline_report.output,
        history=rules.tune.output.history,
        best_hp=rules.tune.output.best_hp,
        model=rules.tune.output.model,
        cell_coverage=expand(
            rules.merge_labels.output.coverage,
            zip,
            sample=cells["sample_id"],
            donor=cells["donor_id"],
            allow_missing=True,
        ),
    output:
        rules.tune.log.notebook.replace("tune", "model_report"),
    conda:
        "../envs/model.yml"
    params:
        random_state=1,
    log:
        notebook=rules.tune.log.notebook.replace("tune", "model_report"),
    notebook:
        "../scripts/model_report.py.ipynb"


rule calls_report:
    input:
        data=rules.tune.output.predictions,
        bulk=rules.germline_report.output,
    output:
        rules.model_report.log.notebook.replace("model_report", "calls_report"),
    conda:
        "../envs/model.yml"
    log:
        notebook=rules.model_report.log.notebook.replace("model_report", "calls_report"),
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
        regions=rules.tune.output.predictions.replace(
            "predictions.pqt", "{donor}.calls.bed"
        ),
        bams=get_donor_bams,
        rmsk=rules.run_rmsk.output.bed,
        megane_gaussian=expand(
            rules.vcf2bed.output, vcf="megane_gaussian", allow_missing=True
        ),
        primer_sites=rules.blast_primers.output.bed,
        fasta=config["genome"]["fasta"],
        fai=config["genome"]["fai"],
    output:
        rules.tune.output.predictions.replace(
            "predictions.pqt", "{donor}.igv_snapshots.bat"
        ),
    log:
        rules.tune.output.predictions.replace(
            "predictions.pqt", "{donor}.igv_snapshots.log"
        ),
    conda:
        "../envs/model.yml"
    params:
        maxPanelHeight=200,
        colorBy="READ_STRAND",
    script:
        "../scripts/igv_snapshots.py"


rule igv_download:
    output:
        "resources/IGV_Linux_2.17.4/igv_hidpi.sh",
    params:
        runtime="600",
        memory="1G",
    shell:
        """
        cd resources
        curl https://data.broadinstitute.org/igv/projects/downloads/2.17/IGV_Linux_2.17.4_WithJava.zip > IGV.zip 2> /dev/null
        unzip -o IGV.zip 2>&1 > /dev/null
        """


rule igv_snapshots:
    input:
        script=rules.make_igv_batch_script.output,
        igv=rules.igv_download.output,
    output:
        directory("{outdir}/results/{genome}/{reads}/peaks/{donor}/snapshots"),
    log:
        "{outdir}/results/{genome}/{reads}/peaks/{donor}/igv_snapshots.log",
    shell:
        """
        # requires xvfb-run to be installed by root for headless display
        xvfb-run -s "-screen 0, 7680x4320x24" --auto-servernum {input.igv} -b {input.script} &> {log}
        """
