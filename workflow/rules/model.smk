rule get_features:
    input:
        bgz=rules.tabix.output.bgz,
        tbi=rules.tabix.output.tbi,
        fa=rules.gen_ref.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        **config["get_features"],
    output:
        f"{outdir}/results/get_features/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.pickle.gz",
    log:
        f"{outdir}/results/get_features/{{ref}}/{{donor}}/{{dna_type}}/{{sample}}.log",
    cache: True
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_features.py"


# define variables relevant to these rules
folds = range(1, config["num_folds"] + 1)
model_ids = list(config["models"].keys())


def get_non_ref_l1(wildcards):
    KNRGL_build = get_KNRGL_build(wildcards)
    if wildcards.ref == "hs37d5":
        return f"{outdir}/resources/{wildcards.db}/{wildcards.ref}_fixnames_insertions.bed"
    elif wildcards.ref != KNRGL_build:
        return f"{outdir}/resources/{wildcards.db}/{wildcards.ref}_lifted_insertions.bed"
    else:
        return f"{outdir}/resources/{wildcards.db}/{wildcards.ref}_insertions.bed"


rule folds:
    input:
        samples=lambda wildcards: expand(
            "{outdir}/results/get_features/{{ref}}/{donor}/{{dna_type}}/{sample}.pickle.gz",
            zip,
            outdir=outdir,
            donor=samples.loc[(samples["dna_type"] == wildcards.dna_type)]["donor"],
            sample=samples.loc[(samples["dna_type"] == wildcards.dna_type)]["sample"],
        ),
        non_ref_l1=get_non_ref_l1,
        ref_l1=rules.run_rmsk.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        num_folds=config["num_folds"],
        min_reads=config["get_features"]["min_reads"],
    output:
        train_features=expand(
            "{outdir}/results/folds/{{ref}}_{{db}}/{{dna_type}}/fold_{fold}/X_train.pickle.gz",
            outdir=outdir,
            fold=folds,
        ),
        test_features=expand(
            "{outdir}/results/folds/{{ref}}_{{db}}/{{dna_type}}/fold_{fold}/X_test.pickle.gz",
            outdir=outdir,
            fold=folds,
        ),
        train_labels=expand(
            "{outdir}/results/folds/{{ref}}_{{db}}/{{dna_type}}/fold_{fold}/Y_train.pickle",
            outdir=outdir,
            fold=folds,
        ),
        test_labels=expand(
            "{outdir}/results/folds/{{ref}}_{{db}}/{{dna_type}}/fold_{fold}/Y_test.pickle",
            outdir=outdir,
            fold=folds,
        ),
        label_encoder=f"{outdir}/results/folds/{{ref}}_{{db}}/{{dna_type}}/label_encoder.pickle",
    log:
        f"{outdir}/results/folds/{{ref}}_{{db}}/{{dna_type}}.log",
    conda:
        "../envs/model.yml"
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
            "{outdir}/results/train_test/{{ref}}_{{db}}/{{dna_type}}/{{model_id}}/fold_{fold}/train_predictions.pickle",
            outdir=outdir,
            fold=folds,
        ),
        train_proba=expand(
            "{outdir}/results/train_test/{{ref}}_{{db}}/{{dna_type}}/{{model_id}}/fold_{fold}/train_probabilities.pickle",
            outdir=outdir,
            fold=folds,
        ),
        test_pred=expand(
            "{outdir}/results/train_test/{{ref}}_{{db}}/{{dna_type}}/{{model_id}}/fold_{fold}/test_predictions.pickle",
            outdir=outdir,
            fold=folds,
        ),
        test_proba=expand(
            "{outdir}/results/train_test/{{ref}}_{{db}}/{{dna_type}}/{{model_id}}/fold_{fold}/test_probabilities.pickle",
            outdir=outdir,
            fold=folds,
        ),
    log:
        f"{outdir}/results/train_test/{{ref}}_{{db}}/{{dna_type}}/{{model_id}}.log",
    conda:
        "../envs/model.yml"
    script:
        "../scripts/train_test.py"


rule metrics:
    input:
        train_labels=rules.folds.output.train_labels,
        train_pred=expand(
            "{outdir}/results/train_test/{{ref}}_{{db}}/{{dna_type}}/{model_id}/fold_{fold}/train_predictions.pickle",
            outdir=outdir,
            model_id=model_ids,
            fold=folds,
        ),
        train_proba=expand(
            "{outdir}/results/train_test/{{ref}}_{{db}}/{{dna_type}}/{model_id}/fold_{fold}/train_probabilities.pickle",
            outdir=outdir,
            model_id=model_ids,
            fold=folds,
        ),
        test_labels=rules.folds.output.test_labels,
        test_pred=expand(
            "{outdir}/results/train_test/{{ref}}_{{db}}/{{dna_type}}/{model_id}/fold_{fold}/test_predictions.pickle",
            outdir=outdir,
            model_id=model_ids,
            fold=folds,
        ),
        test_proba=expand(
            "{outdir}/results/train_test/{{ref}}_{{db}}/{{dna_type}}/{model_id}/fold_{fold}/test_probabilities.pickle",
            outdir=outdir,
            model_id=model_ids,
            fold=folds,
        ),
        label_encoder=rules.folds.output.label_encoder,
    params:
        num_folds=config["num_folds"],
        models=model_ids,
    output:
        prcurve=f"{outdir}/results/metrics/{{ref}}_{{db}}/{{dna_type}}/prcurve.svg",
    conda:
        "../envs/model.yml"
    log:
       f"{outdir}/results/metrics/{{ref}}_{{db}}/{{dna_type}}.log",
    script:
        "../scripts/metrics.py"
