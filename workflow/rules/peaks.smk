if config["genome"]["region"] != "all":
    fai = expand(
        rules.gen_ref.output[1], ref=config["genome"]["build"], outdir=config["outdir"]
    )[0]
    with open(fai, "r") as f:
        genome_size = sum([int(x.split("\t")[1]) for x in f.readlines()])


rule macs2:
    input:
        rules.bwa_mem.output,
    output:
        peaks="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_peaks.narrowPeak",
        summits="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_summits.bed",
        xls="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_peaks.xls",
        cutoff="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_cutoff_analysis.txt",
        bgz="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.bed.gz",
        tbi="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.bed.gz.tbi",
    log:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/peaks.yml"
    params:
        genome_size="hs" if config["genome"]["region"] == "all" else genome_size,
        keep_dup=1,
        qValue_cutoff=0.05
    shell:
        """
        touch {log} && exec > {log} 2>&1

        # get filtered reads
        macs2 filterdup --keep-dup {params.keep_dup} -i {input} | \
            awk '{{if ($2 >= 0 && $3 >= 0) print $0}}' | \
            bedtools sort | \
            bgzip > {output.bgz}
        tabix -p bed {output.bgz}

        macs2 callpeak \
            -g {params.genome_size} \
            -t {input} \
            -q {params.qValue_cutoff} \
            --SPMR \
            --format BAMPE \
            --keep-dup {params.keep_dup} \
            --name {wildcards.sample} \
            --outdir $(dirname {output.peaks}) \
            --cutoff-analysis
        """


rule macs2_evaluate:
    input:
        peaks=rules.macs2.output.peaks,
        reads=rules.macs2.output.bgz,
        ref_l1=rules.run_rmsk.output[0],
    output:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.ipynb",
    log:
        notebook="{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.ipynb",
    conda:
        "../envs/jupyter_peaks.yml"
    notebook:
        "../notebooks/evaluate_macs2.py.ipynb"


rule render_macs2_evaluate:
    input:
        rules.macs2_evaluate.output,
    output:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.html",
    conda:
        "../envs/jupyter_peaks.yml"
    shell:
        "jupyter nbconvert --to html --execute {input} --output $(basename {output})"


rule peaks:
    input:
        expand(
            expand(
                rules.render_macs2_evaluate.output,
                zip,
                dna_type="bulk",
                donor=samples.loc[samples["dna_type"] == "bulk"]["donor"],
                sample=samples.loc[samples["dna_type"] == "bulk"]["sample"],
                allow_missing=True,
            ),
            outdir=config["outdir"],
            ref=config["genome"]["build"],
        ),
