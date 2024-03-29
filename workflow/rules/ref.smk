rule get_line1_consensus:
    output:
        hmm="resources/LINE1_3end/LINE1_3end.hmm",
        fa="resources/LINE1_3end/LINE1_3end.fa",
    log:
        "resources/LINE1_3end/dfam_query.log",
    conda:
        "../envs/ref.yml"
    params:
        accessions=[
            "DF000000225",
            "DF000000339",
            "DF000000340",
            "DF000000341",
            "DF000000342",
            "DF000000343",
            "DF000000344",
            "DF000000345",
            "DF000000347",
            "DF000000327",
            "DF000000329",
            "DF000000331",
            "DF000000333",
            "DF000000335",
            "DF000000336",
            "DF000000337",
        ],
    shell:
        """
        touch {output} {log}
        for a in {params.accessions}; do
            curl -s https://dfam.org/api/families/$a/hmm?format=hmm >> {output.hmm} 2>> {log}
            curl -s https://dfam.org/api/families/$a | jq -r '.name' | awk '{{print ">"$1}}' >> {output.fa} 2>> {log}
            curl -s https://dfam.org/api/families/$a | jq -r '.consensus_sequence' >> {output.fa} 2>> {log}
        done
        """


# generate high accuracy annotation of L1 sequences in references
rule run_rmsk:
    input:
        fa=config["genome"]["fasta"],
        lib=rules.get_line1_consensus.output.hmm,
    output:
        out=multiext(
            config["genome"]["fasta"],
            ".out",
            ".masked",
        ),
        bed=config["genome"]["fasta"] + ".rmsk.bed",
    log:
        config["genome"]["fasta"] + ".rmsk.log",
    conda:
        "../envs/ref.yml"
    params:
        # -s Slow search; 0-5% more sensitive, 2-3 times slower than default;
        # empty string is default
        # -q Quick search; 5-10% less sensitive, 2-5 times faster than default
        # -qq Rush job; about 10% less sensitive, 4->10 times faster than default
        speed="-qq" if config["genome"]["name"] == "chr2122" else "-s",
    shadow:
        "shallow"
    threads: 24
    shell:
        """
        RepeatMasker -pa {threads} -lib {input.lib} -no_is -e hmmer {params.speed} {input} > {log} 2>&1
        rmsk2bed < {output.out[0]} | grep _3end | cut -f 1-15 > {output.bed}
        """


rule blast_primers:
    input:
        ref_fa=config["genome"]["fasta"],
        primer_fa="resources/L1_capture/L1_capture.fa",
    output:
        db=multiext(
            "resources/{genome}/BLASTDB/{genome}",
            ".ndb",
            ".nhr",
            ".nin",
            ".njs",
            ".not",
            ".nsq",
            ".ntf",
            ".nto",
        ),
        bed="resources/{genome}/blast_primers/L1_primers.bed",
    log:
        "resources/{genome}/blast_primers/blast_primers.log",
    conda:
        "../envs/blast.yml"
    script:
        "../scripts/blast_primers.py"


rule en_motif:
    input:
        config["genome"]["fasta"],
    output:
        pos="resources/{genome}/en_pos_score.wig",
        neg="resources/{genome}/en_neg_score.wig",
    log:
        "resources/{genome}/en_motif.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/en_motif.py"
