rule features:
    input:
        bgz=rules.tabix.output.bgz,
        tbi=rules.tabix.output.tbi,
        fa=rules.gen_ref.output[0],
        chromsizes=rules.gen_ref.output[2],
    output:
        multiext(
            "results/features/{ref}/{donor}/{dna_type}/{sample}",
            ".bgz",
            ".bgz.tbi",
            ".header.txt",
            ".unsorted.txt",
            ".sorted.txt",
            ".pickle.gz"
        ),
    log:
        "results/features/{ref}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        touch {log} && exec 2>{log} 

        # get features for each window
        workflow/scripts/get_window_features_occupied.py \
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
            --write_header_to {output[2]} \
            > {output[3]} 

        # Maybe need to add tmpdir?
        # sort features table
        sort --buffer-size=1G -k1,1 -k2,2n -k3,3n < {output[3]} > {output[4]}

        # add header to features table and compress
        cat {output[2]} {output[4]} | bgzip -c > {output[0]}

        tabix -S 1 -s 1 -b 2 -e 3 -0 {output[0]}

        # get additional features: flanking reads
        workflow/scripts/compute_features_and_pickle.py \
            --bgz {output[0]} \
            --chromsizes {input.chromsizes} \
            --outfile {output[5]}
        """


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
