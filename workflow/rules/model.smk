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
        bulk="{outdir}/results/model/get_labels/{ref}_{db}/{donor}_bulk.bed",
        mda="{outdir}/results/model/get_labels/{ref}_{db}/{donor}_mda.pqt",
    log:
        "{outdir}/results/model/get_labels/{ref}_{db}/{donor}.log",
    conda:
        "../envs/features.yml"
    threads: 8
    script:
        "../scripts/get_labels.py"


rule feature_report:
    input:
        expand(
            rules.get_labels.output.mda,
            donor=set(samples["donor"]),
            allow_missing=True,
        ),
    output:
        "{outdir}/results/model/get_labels/{ref}_{db}/feature_report.ipynb",
    log:
        notebook="{outdir}/results/model/get_labels/{ref}_{db}/feature_report.ipynb",
    conda:
        "../envs/model.yml"
    notebook:
        "../notebooks/feature_report.py.ipynb"


rule folds:
    input:
        expand(
            rules.get_labels.output.mda,
            donor=set(samples["donor"]),
            allow_missing=True,
        ),
    params:
        **config["folds"],
    output:
        "{outdir}/results/model/folds/{ref}_{db}/folds.pkl.gz",
    log:
        "{outdir}/results/model/folds/{ref}_{db}/folds.log",
    conda:
        "../envs/model.yml"
    script:
        "../scripts/folds.py"


rule train_test:
    input:
        rules.folds.output,
    params:
        model_name=lambda wc: config["models"][wc.model_id]["name"],
        model_params=lambda wc: config["models"][wc.model_id]["params"],
        train_sampling_strategy=lambda wc: config["models"][wc.model_id][
            "train_sampling_strategy"
        ],
    output:
        "{outdir}/results/model/train_test/{ref}_{db}/{model_id}.pkl.gz",
    log:
        "{outdir}/results/model/train_test/{ref}_{db}/{model_id}.log",
    threads: 8
    benchmark:
        "{outdir}/results/model/train_test/{ref}_{db}/{model_id}.benchmark.txt"
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
        "{outdir}/results/model/train_test/{ref}_{db}/model_report.ipynb",
    conda:
        "../envs/model.yml"
    log:
        notebook="{outdir}/results/model/train_test/{ref}_{db}/model_report.ipynb",
    notebook:
        "../notebooks/model_report.py.ipynb"


rule render_reports:
    input:
        features=rules.feature_report.output,
        model=rules.model_report.output,
    output:
        features="{outdir}/results/model/get_labels/{ref}_{db}/feature_report.html",
        model="{outdir}/results/model/train_test/{ref}_{db}/model_report.html",
    conda:
        "../envs/model.yml"
    log:
        "{outdir}/results/model/train_test/{ref}_{db}/render_reports.log",
    shell:
        """
        touch {log} && exec > {log} 2>&1
        jupyter nbconvert --to html --execute {input.features} --output $(basename {output.features}) 
        jupyter nbconvert --to html --execute {input.model} --output $(basename {output.model}) 
        """
