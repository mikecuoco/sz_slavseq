rule get_features:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        fa=rules.gen_ref.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        **config["get_features"],
    output:
        "{outdir}/results/get_features/{ref}_{db}/{donor}/{dna_type}/{sample}.pickle.gz",
    log:
        "{outdir}/results/get_features/{ref}_{db}/{donor}/{dna_type}/{sample}.log",
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
        "{outdir}/results/get_labels/{ref}_{db}/{donor}.pickle.gz",
    log:
        "{outdir}/results/get_labels/{ref}_{db}/{donor}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_labels.py"


# define variables relevant to these rules
folds = range(1, config["num_folds"] + 1)
model_ids = list(config["models"].keys())


rule folds:
    input:
        samples=expand(
            "{{outdir}}/results/get_labels/{{ref}}_{{db}}/{donor}.pickle.gz",
            donor=set(samples["donor"]),
        ),
    params:
        num_folds=config["num_folds"],
        min_reads=config["get_features"]["min_reads"],
    output:
        train_features=expand(
            "{{outdir}}/results/folds/{{ref}}_{{db}}/fold_{fold}/X_train.pickle.gz",
            fold=folds,
        ),
        test_features=expand(
            "{{outdir}}/results/folds/{{ref}}_{{db}}/fold_{fold}/X_test.pickle.gz",
            fold=folds,
        ),
        train_labels=expand(
            "{{outdir}}/results/folds/{{ref}}_{{db}}/fold_{fold}/Y_train.pickle",
            fold=folds,
        ),
        test_labels=expand(
            "{{outdir}}/results/folds/{{ref}}_{{db}}/fold_{fold}/Y_test.pickle",
            fold=folds,
        ),
        label_encoder="{outdir}/results/folds/{ref}_{db}/label_encoder.pickle",
    log:
        "{outdir}/results/folds/{ref}_{db}/folds.log",
    conda:
        "../envs/model.yml"
    wildcard_constraints:
        dna_type="\w+",
    script:
        "../scripts/folds.py"


rule train_test:
    input:
        train_features=rules.folds.output.train_features,
        train_labels=rules.folds.output.train_labels,
        test_features=rules.folds.output.test_features,
        label_encoder=rules.folds.output.label_encoder,
    params:
        num_folds=config["num_folds"],
        model_params=lambda wc: config["models"][wc.model_id]["params"],
        model_name=lambda wc: config["models"][wc.model_id]["name"],
    output:
        # model="results/train_test/{ref}/{dna_type}/{model}/model.pickle",
        train_pred=expand(
            "{{outdir}}/results/train_test/{{ref}}_{{db}}/{{model_id}}/fold_{fold}/train_predictions.pickle",
            fold=folds,
        ),
        train_proba=expand(
            "{{outdir}}/results/train_test/{{ref}}_{{db}}/{{model_id}}/fold_{fold}/train_probabilities.pickle",
            fold=folds,
        ),
        test_pred=expand(
            "{{outdir}}/results/train_test/{{ref}}_{{db}}/{{model_id}}/fold_{fold}/test_predictions.pickle",
            fold=folds,
        ),
        test_proba=expand(
            "{{outdir}}/results/train_test/{{ref}}_{{db}}/{{model_id}}/fold_{fold}/test_probabilities.pickle",
            fold=folds,
        ),
    log:
        "{outdir}/results/train_test/{ref}_{db}/{model_id}.log",
    conda:
        "../envs/model.yml"
    script:
        "../scripts/train_test.py"


rule metrics:
    input:
        train_labels=rules.folds.output.train_labels,
        train_pred=expand(
            "{{outdir}}/results/train_test/{{ref}}_{{db}}/{model_id}/fold_{fold}/train_predictions.pickle",
            model_id=model_ids,
            fold=folds,
        ),
        train_proba=expand(
            "{{outdir}}/results/train_test/{{ref}}_{{db}}/{model_id}/fold_{fold}/train_probabilities.pickle",
            model_id=model_ids,
            fold=folds,
        ),
        test_labels=rules.folds.output.test_labels,
        test_pred=expand(
            "{{outdir}}/results/train_test/{{ref}}_{{db}}/{model_id}/fold_{fold}/test_predictions.pickle",
            model_id=model_ids,
            fold=folds,
        ),
        test_proba=expand(
            "{{outdir}}/results/train_test/{{ref}}_{{db}}/{model_id}/fold_{fold}/test_probabilities.pickle",
            model_id=model_ids,
            fold=folds,
        ),
        label_encoder=rules.folds.output.label_encoder,
    params:
        num_folds=config["num_folds"],
        models=model_ids,
    output:
        prcurve="{outdir}/results/metrics/{ref}_{db}/prcurve.svg",
    conda:
        "../envs/model.yml"
    log:
        "{outdir}/results/metrics/{ref}_{db}.log",
    script:
        "../scripts/metrics.py"
