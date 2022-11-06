rule get_features:
    input:
        bgz=rules.tabix.output.bgz,
        tbi=rules.tabix.output.tbi,
        fa=rules.gen_ref.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        window_size=config["model"]["window_size"],
        window_step=config["model"]["window_step"],
        min_mapq=40,
        min_ya=20,
        max_yg=15,
        min_secondary_mapq=20,
        library_3_or_5=3,
    output:
        "results/get_features/{ref}/{donor}/{dna_type}/{sample}.pickle.gz",
    log:
        "results/get_features/{ref}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_features.py"


rule folds:
    input:
        samples=expand(
            "results/get_features/{{ref}}/{donor}/{{dna_type}}/{sample}.pickle.gz",
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
        num_folds=config["model"]["num_folds"],
        min_reads=config["model"]["min_reads"],
    output:
        train_features=expand(
            "results/folds/{{ref}}/{{dna_type}}/fold_{fold}/X_train.pickle.gz",
            fold=range(1, config["model"]["num_folds"] + 1),
        ),
        test_features=expand(
            "results/folds/{{ref}}/{{dna_type}}/fold_{fold}/X_test.pickle.gz",
            fold=range(1, config["model"]["num_folds"] + 1),
        ),
        train_labels=expand(
            "results/folds/{{ref}}/{{dna_type}}/fold_{fold}/Y_train.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
        ),
        test_labels=expand(
            "results/folds/{{ref}}/{{dna_type}}/fold_{fold}/Y_test.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
        ),
        label_encoder="results/folds/{ref}/{dna_type}/label_encoder.pickle",
    log:
        "results/folds/{ref}/{dna_type}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/folds.py"


rule train_test:
    input:
        train_features=rules.folds.output.train_features,
        train_labels=rules.folds.output.train_labels,
        test_features=rules.folds.output.test_features,
        label_encoder=rules.folds.output.label_encoder,
    params:
        num_folds=config["model"]["num_folds"],
    output:
        # model="results/train_test/{ref}/{dna_type}/{model}/model.pickle",
        train_predictions=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{{model}}/fold_{fold}/train_predictions.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
        ),
        train_probabilities=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{{model}}/fold_{fold}/train_probabilities.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
        ),
        test_predictions=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{{model}}/fold_{fold}/test_predictions.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
        ),
        test_probabilities=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{{model}}/fold_{fold}/test_probabilities.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
        ),
    log:
        "results/train_test/{ref}/{dna_type}/{model}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/train_test.py"


rule metrics:
    input:
        label_encoder=rules.folds.output.label_encoder,
        train_labels=rules.folds.output.train_labels,
        train_predictions=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{model}/fold_{fold}/train_predictions.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
            model=config["model"]["type"],
        ),
        train_probabilities=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{model}/fold_{fold}/train_probabilities.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
            model=config["model"]["type"],
        ),
        test_labels=rules.folds.output.test_labels,
        test_predictions=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{model}/fold_{fold}/test_predictions.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
            model=config["model"]["type"],
        ),
        test_probabilities=expand(
            "results/train_test/{{ref}}/{{dna_type}}/{model}/fold_{fold}/test_probabilities.pickle",
            fold=range(1, config["model"]["num_folds"] + 1),
            model=config["model"]["type"],
        ),
    params:
        num_folds=config["model"]["num_folds"],
    output:
        prcurve="results/metrics/{ref}/{dna_type}/prcurve.png",
    script:
        "../scripts/metrics.py"
