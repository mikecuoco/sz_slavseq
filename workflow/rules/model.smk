rule get_features:
    input:
        bam=rules.sambamba_sort.output[0],
        bai=rules.sambamba_index.output[0],
        fa=rules.gen_ref.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        **config["get_features"],
    output:
        "{outdir}/results/model/get_features/{donor}/{dna_type}/{sample}.pqt",
    log:
        "{outdir}/results/model/get_features/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_features.py"


def get_labels_input(wildcards):
    donor_samples = samples.loc[samples["donor_id"] == wildcards.donor]
    return {
        "bulk": expand(
            rules.get_features.output,
            sample=donor_samples.loc[samples["dna_type"] == "bulk"]["sample_id"],
            dna_type="bulk",
            allow_missing=True,
        ),
        "mda": expand(
            rules.get_features.output,
            sample=donor_samples.loc[samples["dna_type"] == "mda"]["sample_id"],
            dna_type="mda",
            allow_missing=True,
        ),
    }


rule get_labels:
    input:
        unpack(get_labels_input),
        non_ref_l1=rules.get_donor_knrgl.output[0],
        ref_l1=rules.run_rmsk.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        **config["get_features"],
    output:
        bulk="{outdir}/results/model/get_labels/{donor}_bulk.bed",
        mda="{outdir}/results/model/get_labels/{donor}_mda.pqt",
    log:
        "{outdir}/results/model/get_labels/{donor}.log",
    conda:
        "../envs/features.yml"
    threads: 8
    script:
        "../scripts/get_labels.py"


rule feature_report:
    input:
        expand(
            rules.get_labels.output.mda,
            donor=donors["donor_id"],
            allow_missing=True,
        ),
    output:
        "{outdir}/results/model/get_labels/feature_report.ipynb",
    log:
        notebook="{outdir}/results/model/get_labels/feature_report.ipynb",
    conda:
        "../envs/model.yml"
    notebook:
        "../notebooks/feature_report.py.ipynb"


rule folds:
    input:
        expand(
            rules.get_labels.output.mda,
            donor=donors["donor_id"],
            allow_missing=True,
        ),
    params:
        **config["folds"],
    output:
        folds="{outdir}/results/model/folds/folds.pkl.gz",
        features="{outdir}/results/model/folds/features.txt",
    log:
        "{outdir}/results/model/folds/folds.log",
    conda:
        "../envs/model.yml"
    script:
        "../scripts/folds.py"


rule train_test:
    input:
        folds=rules.folds.output.folds,
        features=rules.folds.output.features,
    params:
        model_name=lambda wc: config["models"][wc.model_id]["name"],
        model_params=lambda wc: config["models"][wc.model_id]["params"],
        train_sampling_strategy=lambda wc: config["models"][wc.model_id][
            "train_sampling_strategy"
        ],
    output:
        "{outdir}/results/model/train_test/{model_id}.pkl.gz",
    log:
        "{outdir}/results/model/train_test/{model_id}.log",
    threads: 8
    benchmark:
        "{outdir}/results/model/train_test/{model_id}.benchmark.txt"
    conda:
        "../envs/model.yml"
    script:
        "../scripts/train_test.py"


rule model_report:
    input:
        expand(
            rules.train_test.output,
            model_id=config["models"].keys(),
            allow_missing=True,
        ),
    output:
        "{outdir}/results/model/train_test/model_report.ipynb",
    conda:
        "../envs/model.yml"
    log:
        notebook="{outdir}/results/model/train_test/model_report.ipynb",
    notebook:
        "../notebooks/model_report.py.ipynb"


rule render_reports:
    input:
        features=rules.feature_report.output,
        model=rules.model_report.output,
    output:
        features="{outdir}/results/model/get_labels/feature_report.html",
        model="{outdir}/results/model/train_test/model_report.html",
    conda:
        "../envs/model.yml"
    log:
        "{outdir}/results/model/train_test/render_reports.log",
    shell:
        """
        touch {log} && exec > {log} 2>&1
        jupyter nbconvert --to html --execute {input.features} --output $(basename {output.features})
        jupyter nbconvert --to html --execute {input.model} --output $(basename {output.model})
        """
