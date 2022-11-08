rule get_features:
    input:
        bgz=rules.tabix.output.bgz,
        tbi=rules.tabix.output.tbi,
        fa=rules.gen_ref.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        **config["get_features"],
    output:
        "results/get_features/{ref}/{donor}/{dna_type}/{sample}.pickle.gz",
    log:
        "results/get_features/{ref}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_features.py"


# define variables relevant to these rules
folds = range(1, config["num_folds"] + 1)
model_ids = list(config["models"].keys())


rule folds:
    input:
        samples=expand(
            "results/get_features/{{ref}}/{donor}/{{dna_type}}/{sample}.pickle.gz",
            zip,
            donor=samples.loc[(samples["dna_type"] == "mda")]["donor"],
            sample=samples.loc[(samples["dna_type"] == "mda")]["sample"],
        ),
        non_ref_l1=expand(
            "resources/{{ref}}/{{ref}}_{db}_insertions.bed",
            db=config["non_ref_germline_l1"]["source"],
        ),
        ref_l1=rules.run_rmsk.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        num_folds=config["num_folds"],
        min_reads=config["get_features"]["min_reads"],
    output:
        train_features=expand(
            "results/folds/{{ref}}/{{dna_type}}/fold_{fold}/X_train.pickle.gz",
            fold=folds,
        ),
        test_features=expand(
            "results/folds/{{ref}}/{{dna_type}}/fold_{fold}/X_test.pickle.gz",
            fold=folds,
        ),
        train_labels=expand(
            "results/folds/{{ref}}/{{dna_type}}/fold_{fold}/Y_train.pickle",
            fold=folds,
        ),
        test_labels=expand(
            "results/folds/{{ref}}/{{dna_type}}/fold_{fold}/Y_test.pickle",
            fold=folds,
        ),
        label_encoder="results/folds/{ref}/{dna_type}/label_encoder.pickle",
    log:
        "results/folds/{ref}/{dna_type}.log",
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
            "results/train_test/{{ref}}/{{dna_type}}/{{model_id}}/fold_{fold}/train_predictions.pickle",
            fold=folds,
        ),
        train_proba=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{{model_id}}/fold_{fold}/train_probabilities.pickle",
            fold=folds,
        ),
        test_pred=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{{model_id}}/fold_{fold}/test_predictions.pickle",
            fold=folds,
        ),
        test_proba=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{{model_id}}/fold_{fold}/test_probabilities.pickle",
            fold=folds,
        ),
    log:
        "results/train_test/{ref}/{dna_type}/{model_id}.log",
    conda:
        "../envs/model.yml"
    script:
        "../scripts/train_test.py"


rule metrics:
    input:
        train_labels=rules.folds.output.train_labels,
        train_pred=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{model_id}/fold_{fold}/train_predictions.pickle",
            model_id=model_ids,
            fold=folds,
        ),
        train_proba=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{model_id}/fold_{fold}/train_probabilities.pickle",
            model_id=model_ids,
            fold=folds,
        ),
        test_labels=rules.folds.output.test_labels,
        test_pred=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{model_id}/fold_{fold}/test_predictions.pickle",
            model_id=model_ids,
            fold=folds,
        ),
        test_proba=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{model_id}/fold_{fold}/test_probabilities.pickle",
            model_id=model_ids,
            fold=folds,
        ),
        label_encoder=rules.folds.output.label_encoder,
    params:
        num_folds=config["num_folds"],
        models=model_ids,
    output:
        prcurve="results/metrics/{ref}/{dna_type}/prcurve.svg",
    conda:
        "../envs/model.yml"
    log:
        "results/metrics/{ref}/{dna_type}.log",
    script:
        "../scripts/metrics.py"
