rule features:
    input:
        bgz=rules.tabix.output.bgz,
        tbi=rules.tabix.output.tbi,
        fa=rules.gen_ref.output[0],
        chromsizes=rules.gen_ref.output[2],
    output:
        # TODO: cleanup with multiext
        bgz="results/features/{ref}/{donor}/{dna_type}/{sample}.bgz",
        tbi="results/features/{ref}/{donor}/{dna_type}/{sample}.bgz.tbi",
        header="results/features/{ref}/{donor}/{dna_type}/{sample}.header.txt",
        unsorted="results/features/{ref}/{donor}/{dna_type}/{sample}.unsorted.txt",
        sorted="results/features/{ref}/{donor}/{dna_type}/{sample}.sorted.txt",
    log:
        "results/features/{ref}/{donor}/{dna_type}/{sample}.log",
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
        chromsizes=rules.gen_ref.output[2],
    output:
        "results/flank_features/{ref}/{donor}/{dna_type}/{sample}.pickle.gz",
    log:
        "results/flank_features/{ref}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/compute_features_and_pickle.py"


rule split_train_test:
    input:
        samples=expand(
            "results/flank_features/{{ref}}/{donor}/{{dna_type}}/{sample}.pickle.gz",
            donor=samples.loc[(samples["dna_type"] == "mda")]["sample"],
            sample=samples.loc[(samples["dna_type"] == "mda")]["donor"],
        ),
        non_ref_l1=expand(
            rules.get_non_ref_l1_windows.output,
            ref=config["genome"]["build"],
            db=config["non_ref_germline_l1"]["source"],
        ),
        ref_l1=rules.get_rmsk_windows.output,
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
