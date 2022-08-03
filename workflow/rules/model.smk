rule features:
    input: 
        bgz = rules.tabix.output.bgz,
        fa = expand(rules.fix_names_clean.output.fa,  ref=config["ref"]),
        chromsizes = expand(rules.fix_names_clean.output.chromsizes,  ref=config["ref"])
    output:
        bgz = "results/features/{donor}/{dna_type}/{sample}.bgz",
        tbi = "results/features/{donor}/{dna_type}/{sample}.bgz.tbi",
        header = "results/features/{donor}/{dna_type}/{sample}.header.txt",
        unsorted = "results/features/{donor}/{dna_type}/{sample}.unsorted.txt",
        sorted = "results/features/{donor}/{dna_type}/{sample}.sorted.txt"
    log: "results/features/{donor}/{dna_type}/{sample}.log"
    conda: "../envs/env.yml"
    shell:
        '''
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
        '''

rule flank_features:
    input: 
        bgz = rules.features.output.bgz,
        chromsizes = expand(rules.fix_names_clean.output.chromsizes, ref=config["ref"])
    output: "results/flank_features/{donor}/{dna_type}/{sample}.pickle.gz"
    log: "results/flank_features/{donor}/{dna_type}/{sample}.log"
    conda: "../envs/env.yml"
    script: "../scripts/compute_features_and_pickle.py"

rule folds:
    input: 
        # TODO: add support to split by type
        samples = expand("results/flank_features/{{donor}}/{dna_type}/{sample}.pickle.gz", 
                        zip,
                        sample=samples['sample'],
                        dna_type=samples['dna_type']), 
        chromsizes = expand(rules.fix_names_clean.output.chromsizes, ref=config["ref"]),
        non_ref_l1 = rules.get_eul1db.output,
        ref_l1 = rules.get_rmsk.output.ref_l1
    params:
        num_folds = config["model"]["num_folds"],
        min_reads = config["model"]["min_reads"],
        window_size = config["model"]["window_size"]
    output: multiext("results/folds/{donor}/", "X_train.pickle.gz", "X_test.pickle.gz", "Y_train.pickle.gz", "Y_test.pickle.gz")
    log: "results/folds/{donor}/log.log"
    conda: "../envs/env.yml"
    script: "../scripts/folds.py"