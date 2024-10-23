rule make_bigwig:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        chrom_sizes=config["genome"]["genome"],
    output:
        bg=temp("{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.bg"),
        bw="{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.bw",
    log:
        "{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.make_bigwig.log",
    conda:
        "../envs/ref.lock.yml"
    shell:
        """
        # for reads that are: not read2, not secondary, not duplicate, not unmapped, generated bedgraph of 3' ends
        samtools view -h -F 1412 {input.bam} | bedtools genomecov -bg -ibam stdin -3 | sort -k1,1 -k2,2n > {output.bg} 2> {log}
        # convert bedgraph to bigwig
        bedGraphToBigWig {output.bg} {input.chrom_sizes} {output.bw} 2>> {log}
        """


rule merge_bigwig_donors:
    input:
        bw=lambda wc: expand(
            rules.make_bigwig.output.bw,
            sample=cells.loc[cells["donor_id"] == wc.donor, "sample_id"],
            allow_missing=True,
        ),
        chrom_sizes=config["genome"]["genome"],
    output:
        bg=temp("{outdir}/results/{genome}/{reads}/bigwig/{donor}.bg"),
        bg_sorted=temp("{outdir}/results/{genome}/{reads}/bigwig/{donor}.sorted.bg"),
        bw="{outdir}/results/{genome}/{reads}/bigwig/{donor}.bw",
    log:
        "{outdir}/results/{genome}/{reads}/bigwig/{donor}.merge_bigwig.log",
    conda:
        "../envs/ref.lock.yml"
    shell:
        """
        bigWigMerge -max -threshold=10 {input.bw} {output.bg} 2> {log}
        sort -k1,1 -k2,2bn {output.bg} > {output.bg_sorted} 2>> {log}
        bedGraphToBigWig {output.bg_sorted} {input.chrom_sizes} {output.bw} 2>> {log}
        """


def get_group(wildcards):
    if wildcards.group == "bulk":
        return {
            "sample": bulk["sample_id"],
            "donor": bulk["donor_id"],
        }
    elif wildcards.group == "cells":
        return {
            "sample": cells["sample_id"],
            "donor": cells["donor_id"],
        }
    else:
        raise ValueError(f"Unknown sample group: {wildcards.group}")


rule merge_bigwig_all:
    input:
        bw=lambda wc: expand(
            rules.make_bigwig.output.bw,
            zip,
            **get_group(wc),
            allow_missing=True,
        ),
        chrom_sizes=config["genome"]["genome"],
    output:
        bg=temp("{outdir}/results/{genome}/{reads}/bigwig/{group}.bg"),
        bg_sorted=temp("{outdir}/results/{genome}/{reads}/bigwig/{group}.sorted.bg"),
        bw="{outdir}/results/{genome}/{reads}/bigwig/{group}.bw",
    log:
        "{outdir}/results/{genome}/{reads}/bigwig/{group}.merge_bigwig.log",
    conda:
        "../envs/ref.lock.yml"
    shell:
        """
        bigWigMerge -max -threshold=10 {input.bw} {output.bg} 2> {log}
        sort -k1,1 -k2,2bn {output.bg} > {output.bg_sorted} 2>> {log}
        bedGraphToBigWig {output.bg_sorted} {input.chrom_sizes} {output.bw} 2>> {log}
        """


rule make_bw_peaks:
    input:
        bw="{outdir}/results/{genome}/{reads}/bigwig/{donor}/{sample}.bw",
    output:
        bed="{outdir}/results/{genome}/{reads}/bw_peaks/{donor}/{sample}.peaks.bed",
    log:
        "{outdir}/results/{genome}/{reads}/bw_peaks/{donor}/{sample}.make_peaks.log",
    conda:
        "../envs/features.yml"
    params:
        size=200,
        step=1,
        minreads=5,
    script:
        "../scripts/make_bw_peaks.py"


rule make_bw_cons_peaks:
    input:
        bw=rules.merge_bigwig_all.output.bw,
    output:
        bed="{outdir}/results/{genome}/{reads}/bw_cons_peaks/{group}.peaks.bed",
    log:
        "{outdir}/results/{genome}/{reads}/bw_cons_peaks/{group}.make_peaks.log",
    conda:
        "../envs/features.yml"
    params:
        size=200,
        step=1,
        minreads=10,
    script:
        "../scripts/make_bw_peaks.py"


# rule make_bam_peaks:
# 	input:
# 		bam=rules.sort.output[0],
# 		bai=rules.index.output[0],
# 		chrom_sizes=config["genome"]["genome"],
# 	output:
# 		bed="{outdir}/results/{genome}/{reads}/peaks/{donor}/{sample}.bed.gz",
# 		tbi="{outdir}/results/{genome}/{reads}/peaks/{donor}/{sample}.bed.gz.tbi",
# 	log:
# 		"{outdir}/results/{genome}/{reads}/peaks/{donor}/{sample}.make_bam_peaks.log",
# 	conda:
# 		"../envs/features.yml"
# 	script:
# 		"../scripts/make_bam_peaks.py"

# rule make_breakpoints:
# 	input:
# 		bam=rules.sort.output[0],
# 		bai=rules.index.output[0],
# 		chrom_sizes=config["genome"]["genome"],
# 	output:
# 		bed="{outdir}/results/{genome}/{reads}/breakpoints/{donor}/{sample}.bed.gz",
# 		tbi="{outdir}/results/{genome}/{reads}/breakpoints/{donor}/{sample}.bed.gz.tbi",
# 	log:
# 		"{outdir}/results/{genome}/{reads}/breakpoints/{donor}/{sample}.make_breakpoints.log",
# 	conda:
# 		"../envs/features.yml"
# 	script:
# 		"../scripts/make_breakpoints.py"


rule make_greedy_peaks:
    input:
        bw=rules.merge_bigwig_donors.output.bw,
    output:
        bed="{outdir}/results/{genome}/{reads}/greedy/{donor}.peaks.bed",
    log:
        "{outdir}/results/{genome}/{reads}/greedy/{donor}_make_greedy_peaks.log",
    conda:
        "../envs/features.yml"
    params:
        blacklist_bp=40000,
        extend_bp=250,
    script:
        "../scripts/make_greedy_peaks.py"


def get_peaks(wildcards):
    if wildcards.region == "greedy":
        return rules.make_greedy_peaks.output.bed
    elif wildcards.region == "bw_peaks":
        return rules.make_bw_peaks.output.bed
    elif wildcards.region == "bw_cons_peaks":
        return rules.make_bw_cons_peaks.output.bed
    else:
        raise ValueError(f"Unknown region: {wildcards.region}")


rule get_peak_features:
    input:
        bed=get_peaks,
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
    output:
        bed="{outdir}/results/{genome}/{reads}/{region}/{donor}/{sample}.features.bed.gz",
        tbi="{outdir}/results/{genome}/{reads}/{region}/{donor}/{sample}.features.bed.gz.tbi",
    log:
        "{outdir}/results/{genome}/{reads}/{region}/{donor}/{sample}.get_features.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_peak_features.py"


rule merge_regions:
    input:
        lambda wc: expand(
            rules.get_peak_features.output.bed,
            zip,
            **get_group(wc),
            allow_missing=True,
        ),
    output:
        bed="{outdir}/results/{genome}/{reads}/{region}/{group}.features.bed.gz",
        tbi="{outdir}/results/{genome}/{reads}/{region}/{group}.features.bed.gz.tbi",
    log:
        "{outdir}/results/{genome}/{reads}/{region}/{group}_merge.log",
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


labels = [
    "bulk",
    "megane",
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
]


# TODO: add closest and overlap labels
rule label_regions:
    input:
        **rmsk_anno,
        **knrgl_anno,
        bulk="{outdir}/results/{genome}/{reads}/{region}/bulk.features.bed.gz",
        cell_regions="{outdir}/results/{genome}/{reads}/{region}/{group}.features.bed.gz",
        primer_sites=rules.blast_primers.output.bed,
    output:
        bed="{outdir}/results/{genome}/{reads}/{region}/{group}.labelled.bed.gz",
        tbi="{outdir}/results/{genome}/{reads}/{region}/{group}.labelled.bed.gz.tbi",
        pqt="{outdir}/results/{genome}/{reads}/{region}/{group}.labelled.pqt",
    log:
        rules.merge_regions.log[0].replace("merge", "label"),
    run:
        from subprocess import run
        from pathlib import Path
        import pandas as pd

        # make log file
        if Path(str(log)).exists():
            Path(str(log)).unlink()
        run(f"touch {log}", shell=True)

        # get header
        header = (
            run(
                f"zcat {input.cell_regions} | head -n 1",
                shell=True,
                capture_output=True,
            )
            .stdout.decode()
            .rstrip("\n")
        )
        cmd = f"bedtools annotate -counts -i {input.cell_regions} -files "
        for l in labels:
            cmd += str(input[l]) + " "
            header += f"\t{l}"

            # make temp unzipped bed output, add header
        bed = str(output.bed).replace(".gz", "")
        run(f"echo '{header}' > {bed} 2>> {log}", shell=True)

        # add data
        cmd += f"| sort -k1,1 -k2,2n >> {bed} 2>> {log}"
        run(cmd, shell=True)

        # compress
        run(f"bgzip {bed} && tabix -f -p bed {output.bed} 2>> {log}", shell=True)

        # read it into mem to save to more efficient formats
        data = pd.read_csv(output.bed, sep="\t")
        data.columns = data.columns.str.replace("#", "")
        data["locus"] = (
            data["Chromosome"].astype(str)
            + ":"
            + data["Start"].astype(str)
            + "-"
            + data["End"].astype(str)
        )
        data["donor_id"] = data["donor_id"].astype(str)
        data["width"] = data["End"] - data["Start"]
        data["tissue"] = data["cell_id"].apply(
            lambda x: "DLPFC" if "usd" in x.lower() else "HIP"
        )
        data[labels] = data[labels].astype(bool)

        # save to pqt
        data.to_parquet(output.pqt)


rule bulk_regions_report:
    input:
        bulk=expand(
            rules.label_regions.output.bed,
            region=["bw_peaks", "greedy", "bw_cons_peaks"],
            group="bulk",
            allow_missing=True,
        ),
        l1hs_rmsk=rmsk_anno["l1hs"],
        megane=expand(
            rules.merge_bed.output,
            vcf="megane",
            allow_missing=True,
        ),
    output:
        "{outdir}/results/{genome}/{reads}/bulk.ipynb",
    log:
        notebook="{outdir}/results/{genome}/{reads}/bulk.ipynb",
    conda:
        "../envs/model.yml"
    notebook:
        "../scripts/bulk_regions_report.py.ipynb"


rule wrangle_regions:
    input:
        **knrgl_anno,
        pqt=expand(rules.label_regions.output.pqt, group="cells", allow_missing=True),
        l1hs_rmsk=rmsk_anno["l1hs"],
        fasta=config["genome"]["fasta"],
        bulk=expand(rules.label_regions.output.bed, group="bulk", allow_missing=True),
        flasch="resources/NIHMS1523125-supplement-8.xls",
    output:
        "{outdir}/results/{genome}/{reads}/{region}/data_ready.pkl",
    log:
        "{outdir}/results/{genome}/{reads}/{region}/wrangle_regions.log",
    conda:
        "../envs/model.yml"
    threads: 12  # dont make this too high or will run out of memory
    script:
        "../scripts/wrangle_regions.py"


rule cell_regions_report:
    input:
        regions=expand(
            rules.wrangle_regions.output,
            region=["greedy", "bw_peaks", "bw_cons_peaks"],
            allow_missing=True,
        ),
        apua="resources/meiyan_apua_validation/U01_L1_primers.chm13v2.csv",
    output:
        "{outdir}/results/{genome}/{reads}/calls.ipynb",
    conda:
        "../envs/model.yml"
    log:
        notebook="{outdir}/results/{genome}/{reads}/calls.ipynb",
    notebook:
        "../scripts/regions_report.py.ipynb"


rule model:
    input:
        data="{outdir}/results/{genome}/{reads}/{region}/data.bed.gz",
    output:
        model="{outdir}/results/{genome}/{reads}/{region}/model.pkl",
        best_hp="{outdir}/results/{genome}/{reads}/{region}/best_hp.json",
        history="{outdir}/results/{genome}/{reads}/{region}/history.history",
        predictions="{outdir}/results/{genome}/{reads}/{region}/predictions.pqt",
        shap_values="{outdir}/results/{genome}/{reads}/{region}/shap_values.pqt",
    conda:
        "../envs/model_gpu.yml" if shutil.which("nvidia-smi") else "../envs/model.yml"
    params:
        max_iter=50,
        random_state=1,
        # features=lambda wc: feature_sets[wc.features],
        # TODO: add list of models to test
        features=[
            "n_reads",
            "n_duplicates",
            "n_unique_5end",
            "three_end_clippedA_mean",
            "germline_distance",
        ],
    threads: 1
    log:
        notebook="{outdir}/results/{genome}/{reads}/{region}/tune.ipynb",
    notebook:
        "../scripts/tune.py.ipynb"


# rule calls_report:
#     input:
#         data=rules.model.output.predictions,
#     output:
#         "{outdir}/results/{genome}/{reads}/{region}/calls.ipynb",
#     conda:
#         "../envs/model.yml"
#     log:
#         notebook="{outdir}/results/{genome}/{reads}/{region}/calls.ipynb",
#     notebook:
#         "../scripts/calls_report.py.ipynb"


rule regions:
    input:
        expand(
            rules.cell_regions_report.output,
            reads="filtered",
            genome=config["genome"]["name"],
            outdir=config["outdir"],
        ),
