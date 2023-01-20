rule macs2:
    input:
        rules.tags.output,
    output:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}_peaks.narrowPeak",
    log:
        "{outdir}/results/macs2/{ref}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/peaks.yml"
    shell:
        "macs2 callpeak -t {input} --name {wildcards.sample} --outdir $(dirname {output}) --nomodel --extsize 750 2> {log}"


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
            rules.tags.output,
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
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/evaluate_macs2.py.ipynb"


rule peaks:
    input:
        expand(
            rules.macs2_evaluate.output,
            outdir=config["outdir"],
            ref=config["genome"]["build"],
            db=list(config["KNRGL"].keys()),
            donor=samples["donor"],
            allow_missing=True,
        ),
