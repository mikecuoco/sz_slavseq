rule features:
    input:
        bgz=rules.tabix.output.bgz,
        tbi=rules.tabix.output.tbi,
        fa=expand(rules.fix_names_clean.output.fa, ref=config["ref"]["build"]),
        chromsizes=expand(
            rules.fix_names_clean.output.chromsizes, ref=config["ref"]["build"]
        ),
    output:
        bgz="results/features/{donor}/{dna_type}/{sample}.bgz",
        tbi="results/features/{donor}/{dna_type}/{sample}.bgz.tbi",
        header="results/features/{donor}/{dna_type}/{sample}.header.txt",
        unsorted="results/features/{donor}/{dna_type}/{sample}.unsorted.txt",
        sorted="results/features/{donor}/{dna_type}/{sample}.sorted.txt",
    log:
        "results/features/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        touch {log} && exec 2>{log} 

        pyslavseq_extract_features \
            --genome_fasta_file {input.fa} \
            --library_3_or_5 3 \
            --occupied \
            --min_mapq 40 \
            --min_ya 20 \
            --max_yg 15 \
            --chromsizes {input.chromsizes} \
            --window_size 750 \
            --window_step 250 \
            --min_secondary_mapq 20 \
            {input.bgz} \
            --write_header_to {output.header} \
            > {output.unsorted} 

        # Maybe need to add tmpdir?
        sort --buffer-size=1G -k1,1 -k2,2n -k3,3n < {output.unsorted} > {output.sorted}

        cat {output.header} {output.sorted} | bgzip -c > {output.bgz}

        tabix -S 1 -s 1 -b 2 -e 3 -0 {output.bgz}
        """


rule flank_features:
    input:
        bgz=rules.features.output.bgz,
        tbi=rules.features.output.tbi,
        chromsizes=expand(
            rules.fix_names_clean.output.chromsizes, ref=config["ref"]["build"]
        ),
    output:
        "results/flank_features/{donor}/{dna_type}/{sample}.pickle.gz",
    log:
        "results/flank_features/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/compute_features_and_pickle.py"


rule folds:
    input:
        samples=get_folds_input_samples,
        chromsizes=expand(
            rules.fix_names_clean.output.chromsizes, ref=config["ref"]["build"]
        ),
        non_ref_l1=get_non_ref_l1,
        ref_l1=expand(rules.get_rmsk.output.ref_l1, ref=config["ref"]["build"]),
    params:
        # non_ref_db=config["ref"]["database"]
        num_folds=config["model"]["num_folds"],
        min_reads=config["model"]["min_reads"],
        fold_window=config["model"]["fold_window"],
    output:
        expand(
            "results/folds/{{donor}}/{{dna_type}}/{fold}/{file}",
            fold=fold_dirs,
            file=[
                "X_train.pickle.gz",
                "X_test.pickle.gz",
                "Y_train.pickle.gz",
                "Y_test.pickle.gz",
            ],
        ),
    log:
        "results/folds/{donor}/{dna_type}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/folds.py"


rule train_test:
    input:
        expand(
            "results/folds/{{donor}}/{{dna_type}}/{fold}/{file}",
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
            "results/train_test/{{donor}}/{{dna_type}}/{fold}/Testing_y_pred.csv",
            fold=fold_dirs,
        ),
    params:
        num_folds=config["model"]["num_folds"],
    log:
        "results/train_test/{donor}/{dna_type}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/rfc.py"


rule somatic_summary:
    input:
        rules.train_test.output,
    output:
        expand(
            "results/somatic_summary/{{donor}}/{{dna_type}}/{file}",
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
        "results/somatic_summary/{donor}/{dna_type}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/somatic_summary.py"
