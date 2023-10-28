rule get_l1hs_consensus:
    output:
        "resources/LINE1/LINE1.fa",
    log:
        "resources/LINE1/dfam_query.log",
    conda:
        "../envs/ref.yml"
    params:
        accessions=["DF0000225"],
    shell:
        """
        touch {log} && exec 1>{log} 2>&1
        touch {output}
        for a in {params.accessions}; do
            curl -s https://dfam.org/api/families/$a | jq -r '.name' | awk '{{print ">"$1}}'  >> {output}
            curl -s https://dfam.org/api/families/$a | jq -r '.consensus_sequence' >> {output}
        done
        """


rule get_line1_hmm:
    output:
        "resources/LINE1_3end/LINE1_3end.hmm",
    log:
        "resources/LINE1_3end/dfam_query.log",
    conda:
        "../envs/ref.yml"
    params:
        accessions=[
            "DF0000225",
            "DF0000339",
            "DF0000340",
            "DF0000341",
            "DF0000342",
            "DF0000343",
            "DF0000344",
            "DF0000345",
            "DF0000347",
            "DF0000327",
            "DF0000329",
            "DF0000331",
            "DF0000333",
            "DF0000335",
            "DF0000336",
            "DF0000337",
        ],
    shell:
        """
        touch {log} && exec 1>{log} 2>&1
        touch {output}
        for a in {params.accessions}; do
            curl -s https://dfam.org/api/families/$a/hmm?format=hmm >> {output}
        done
        """


# generate high accuracy annotation of L1 sequences in references
rule run_rmsk:
    input:
        fa=config["genome"]["fasta"],
        lib=rules.get_line1_hmm.output,
    output:
        multiext(
            config["genome"]["fasta"],
            ".out",
            ".masked",
        ),
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


rule rmsk_to_bed:
    input:
        rules.run_rmsk.output[0],
    output:
        rmsk="{outdir}/resources/rmsk.bed",
        rmsk_1kb_3end="{outdir}/resources/rmsk_1kb_3end.bed",
    log:
        "{outdir}/resources/rmsk_to_bed.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/rmsk_to_bed.py"


# get vcf, convert to bed, remove orphan insertions (higher FP rate and won't be detected in SLAV-seq)
rule xtea_to_bed:
    input:
        lambda wc: donors.loc[wc.donor]["xtea"],
    output:
        xtea="{outdir}/resources/{donor}_insertions.bed",
        xtea_1kb_3end="{outdir}/resources/{donor}_insertions_1kb_3end.bed",
    conda:
        "../envs/features.yml"
    log:
        "{outdir}/resources/{donor}_to_bed.log",
    script:
        "../scripts/xtea_to_bed.py"
