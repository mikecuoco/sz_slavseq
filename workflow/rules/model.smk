rule get_features:
    input:
        bgz=rules.tabix.output.bgz,
        tbi=rules.tabix.output.tbi,
        fa=rules.gen_ref.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        **config["get_features"],
    output:
        "{outdir}/results/model/get_features/{ref}_{db}/{donor}/{dna_type}/{sample}.pqt",
    log:
        "{outdir}/results/model/get_features/{ref}_{db}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_features.py"


def get_non_ref_l1(wildcards):
    KNRGL_build = get_KNRGL_build(wildcards)
    if wildcards.ref == "hs37d5":
        return f"{wildcards.outdir}/resources/{wildcards.db}/{wildcards.ref}_fixnames_insertions.bed"
    elif wildcards.ref != KNRGL_build:
        return f"{wildcards.outdir}/resources/{wildcards.db}/{wildcards.ref}_lifted_insertions.bed"
    else:
        return f"{wildcards.outdir}/resources/{wildcards.db}/{wildcards.ref}_insertions.bed"


def get_labels_input(wildcards):
    donor_samples = samples.loc[samples["donor"] == wildcards.donor]
    return {
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
        "features": expand(
            rules.get_features.output,
            sample=donor_samples.loc[samples["dna_type"] == "mda"]["sample"],
            dna_type="mda",
            allow_missing=True,
        ),
    }


rule get_labels:
    input:
        unpack(get_labels_input),
        non_ref_l1=get_non_ref_l1,
        ref_l1=rules.run_rmsk.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        **config["get_features"],
    output:
        bulk="{outdir}/results/model/get_labels/{ref}_{db}/{label_config}/{donor}_bulk.bed",
        mda="{outdir}/results/model/get_labels/{ref}_{db}/{label_config}/{donor}_mda.pqt",
    log:
        "{outdir}/results/model/get_labels/{ref}_{db}/{label_config}/{donor}.log",
    conda:
        "../envs/features.yml"
    threads: 8
    script:
        "../scripts/get_labels.py"


rule folds:
    input:
        samples=expand(
            rules.get_labels.output.mda,
            donor=set(samples["donor"]),
            allow_missing=True,
        ),
    params:
        **config["folds"],
        min_reads=config["get_features"]["min_reads"],
    output:
        features="{outdir}/results/model/folds/{ref}_{db}/{label_config}/features.pickle",
        labels="{outdir}/results/model/folds/{ref}_{db}/{label_config}/labels.pickle",
        label_encoder="{outdir}/results/model/folds/{ref}_{db}/{label_config}/label_encoder.pickle",
        folds="{outdir}/results/model/folds/{ref}_{db}/{label_config}/folds.pickle.gz",
    log:
        "{outdir}/results/model/folds/{ref}_{db}/{label_config}/folds.log",
    conda:
        "../envs/model.yml"
    threads: 8
    script:
        "../scripts/folds.py"


rule train_test:
    input:
        features=rules.folds.output.features,
        labels=rules.folds.output.labels,
        label_encoder=rules.folds.output.label_encoder,
    params:
        model_params=lambda wc: config["models"][wc.model_id]["params"],
        model_name=lambda wc: config["models"][wc.model_id]["name"],
    output:
        # model="results/train_test/{ref}/{dna_type}/{model}/model.pickle",
        pred="{outdir}/results/model/train_test/{ref}_{db}/{label_config}/{model_id}/pred.pickle",
        proba="{outdir}/results/model/train_test/{ref}_{db}/{label_config}/{model_id}/proba.pickle",
        metrics="{outdir}/results/model/train_test/{ref}_{db}/{label_config}/{model_id}/metrics.pickle",
    log:
        "{outdir}/results/model/train_test/{ref}_{db}/{label_config}/{model_id}/train_test.log",
    threads: 8
    benchmark:
        "{outdir}/results/model/train_test/{ref}_{db}/{label_config}/{model_id}/train_test.benchmark.txt"
    conda:
        "../envs/model.yml"
    script:
        "../scripts/train_test.py"


rule model_report:
    input:
        folds=rules.folds.output.folds,
        metrics=expand(
            rules.train_test.output.metrics,
            model_id=list(config["models"].keys()),
            allow_missing=True,
        ),
        label_encoder=rules.folds.output.label_encoder,
    output:
        "{outdir}/results/model/train_test/{ref}_{db}/{label_config}/model_report.ipynb",
    conda:
        "../envs/jupyter.yml"
    params:
        model_ids=list(config["models"].keys()),
    log:
        notebook="{outdir}/results/model/train_test/{ref}_{db}/{label_config}/model_report.ipynb",
    notebook:
        "../notebooks/model_report.py.ipynb"


rule render_report:
    input:
        notebook=rules.model_report.output,
    output:
        "{outdir}/results/model/train_test/{ref}_{db}/{label_config}/model_report.html",
    conda:
        "../envs/jupyter.yml"
    log:
        "{outdir}/results/model/train_test/{ref}_{db}/{label_config}/model_report.log",
    shell:
        "jupyter nbconvert --to html --execute {input.notebook} --output $(basename {output}) 2> {log} 2>&1"
