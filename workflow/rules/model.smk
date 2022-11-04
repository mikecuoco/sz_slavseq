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


rule split_train_test:
    input:
        samples=expand(
            "results/flank_features/{{ref}}/{donor}/{{dna_type}}/{sample}.pickle.gz",
            donor=samples.loc[(samples["dna_type"] == "mda")]["sample"],
            sample=samples.loc[(samples["dna_type"] == "mda")]["donor"],
        ),
        non_ref_l1=expand(
            "resources/{{ref}}/{{ref}}_{db}_insertions.bed",
            db=config["non_ref_germline_l1"]["source"],
        ),
        ref_l1=rules.run_rmsk.output[0],
    params:
        num_folds=config["model"]["num_folds"],
        min_reads=config["model"]["min_reads"],
        fold_window=config["model"]["fold_window"],
    output:
        train="results/split_train_test/{ref}/{dna_type}/Training_y_pred.csv",
        test="results/split_train_test/{ref}/{dna_type}/Testing_y_pred.csv",
    log:
        "results/split_train_test/{ref}/{dna_type}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/split_train_test.py"
