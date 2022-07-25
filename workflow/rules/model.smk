rule features:
    input: 
        bgz = rules.tabix.output.bgz,
        fa = expand(rules.fix_names_clean.output.fa,  ref=config["ref"]),
        chromsizes = expand(rules.fix_names_clean.output.chromsizes,  ref=config["ref"])
    output:
        bgz = "results/features/{sample}/{donor}_{type}.bgz",
        tbi = "results/features/{sample}/{donor}_{type}.bgz.tbi",
        header = "results/features/{sample}/{donor}_{type}.header.txt",
        unsorted = "results/features/{sample}/{donor}_{type}.unsorted.txt",
        sorted = "results/features/{sample}/{donor}_{type}.sorted.txt"
    log: "results/features/{sample}/{donor}_{type}.log"
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
    output: "results/flank_features/{sample}/{donor}_{type}.pickle.gz"
    log: "results/flank_features/{sample}/{donor}_{type}.log"
    conda: "../envs/env.yml"
    script: "../scripts/compute_features_and_pickle.py"
