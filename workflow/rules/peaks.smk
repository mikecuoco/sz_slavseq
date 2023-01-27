rule macs2:
    input:
        bam=rules.tags.output,
        fai=rules.gen_ref.output[1]
    output:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_peaks.narrowPeak",
    log:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/peaks.yml"
    shell:
        """
        touch {log} && exec > {log} 2>&1
        GSIZE=$(awk '{{sum+=$2}} END{{print sum}}' {input.fai})
        macs2 callpeak -g $GSIZE -t {input.bam} --name {wildcards.sample} --outdir $(dirname {output}) --nomodel --extsize 750 2> {log}
        """


def get_evaluate_input(wildcards):
    donor_samples = samples.loc[samples["donor"] == wildcards.donor]
    return {
        "peaks": expand(
            rules.macs2.output,
            sample=donor_samples.loc[samples["dna_type"] == "bulk"]["sample"],
            dna_type="bulk",
            allow_missing=True,
        ),
        "bam": expand(
            rules.sort.output,
            sample=donor_samples.loc[samples["dna_type"] == "bulk"]["sample"],
            dna_type="bulk",
            allow_missing=True,
        ),
        "bai": expand(
            rules.index.output,
            sample=donor_samples.loc[samples["dna_type"] == "bulk"]["sample"],
            dna_type="bulk",
            allow_missing=True,
        ),
        "bgz": expand(
            rules.tabix.output.bgz,
            sample=donor_samples.loc[samples["dna_type"] == "bulk"]["sample"],
            dna_type="bulk",
            allow_missing=True,
        ),
        "tbi": expand(
            rules.tabix.output.tbi,
            sample=donor_samples.loc[samples["dna_type"] == "bulk"]["sample"],
            dna_type="bulk",
            allow_missing=True,
        ),
    }


rule macs2_evaluate:
    input:
        unpack(get_evaluate_input),
        labels=expand(rules.get_labels.output.bulk, label_config="mergeRL1", allow_missing=True),
        non_ref_l1=get_non_ref_l1,
        ref_l1=rules.run_rmsk.output[0],
    output:
        "{outdir}/results/macs2_eval/{ref}_{db}/{donor}.ipynb",
    log:
        notebook="{outdir}/results/macs2_eval/{ref}_{db}/{donor}.ipynb",
    conda:
        "../envs/jupyter_peaks.yml"
    notebook:
        "../notebooks/evaluate_macs2.py.ipynb"

rule render_macs2_evaluate:
    input:
        rules.macs2_evaluate.output,
    output:
        "{outdir}/results/macs2_eval/{ref}_{db}/{donor}.html",
    conda:
        "../envs/jupyter_peaks.yml"
    shell:
        "jupyter nbconvert --to html --execute {input} --output $(basename {output})"


rule peaks:
    input:
        expand(
            rules.render_macs2_evaluate.output,
            outdir=config["outdir"],
            ref=config["genome"]["build"],
            db=list(config["KNRGL"].keys()),
            donor=samples["donor"],
            allow_missing=True,
        ),
