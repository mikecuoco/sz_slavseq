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
        samples=get_folds_input_samples,
        chromsizes=rules.gen_ref.output[2],
        non_ref_l1=non_ref_l1_bed,
        ref_l1=rules.run_rmsk.output[0],
    params:
        num_folds=config["model"]["num_folds"],
        min_reads=config["model"]["min_reads"],
        fold_window=config["model"]["fold_window"],
    output:
        expand(
            "results/folds/{{ref}}/{{donor}}/{{dna_type}}/{fold}/{file}",
            fold=fold_dirs,
            file=[
                "X_train.pickle.gz",
                "X_test.pickle.gz",
                "Y_train.pickle.gz",
                "Y_test.pickle.gz",
            ],
        ),
    log:
        "results/folds/{ref}/{donor}/{dna_type}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/folds.py"


rule train_test:
    input:
        expand(
            "results/folds/{{ref}}/{{donor}}/{{dna_type}}/{fold}/{file}",
            fold=fold_dirs,
            file=[
                "X_train.pickle.gz",
                "X_test.pickle.gz",
                "Y_train.pickle.gz",
                "Y_test.pickle.gz",
            ],
        ),
    output:
        expand(
            "results/train_test/{{ref}}/{{donor}}/{{dna_type}}/{fold}/Testing_y_pred.csv",
            fold=fold_dirs,
        ),
    params:
        num_folds=config["model"]["num_folds"],
    log:
        "results/train_test/{ref}/{donor}/{dna_type}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/rfc.py"


rule somatic_summary:
    input:
        rules.train_test.output,
    output:
        expand(
            "results/somatic_summary/{{ref}}/{{donor}}/{{dna_type}}/{file}",
            file=[
                "Merged_y_pred.csv",
                "slavseq_sz-intersections-cluster.csv",
                "somatic_candidates-cluster.csv",
                "Cross_validation_metrics.csv",
            ],
        ),
    params:
        num_folds=config["model"]["num_folds"],
        min_reads=config["model"]["min_reads"],
        window_size=config["model"]["window_size"],
        min_prob=config["model"]["prob"],
    log:
        "results/somatic_summary/{ref}/{donor}/{dna_type}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/somatic_summary.py"
