rule install_bwakit:
    output:
        directory("resources/bwa.kit"),
    conda:
        "../envs/env.yml"
    log:
        "resources/install_bwakit.log",
    shell:
        """
        mkdir -p resources && cd resources
        wget -O- -q --no-config https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2 | tar xfj -
        """


rule gen_ref:
    input:
        rules.install_bwakit.output,
    output:
        "resources/{ref}/genome_og.fa",
    log:
        "resources/{ref}/gen_ref.log",
    conda:
        "../envs/env.yml"
    params:
        region=config["genome"]["region"],
    cache: True
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        # run bwa.kit function
        if [ {wildcards.ref} == "hs38DH" ]; then
            wget -q --no-config ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
            gunzip GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
            cat GCA_000001405.15_GRCh38_full_analysis_set.fna {input}/resource-GRCh38/hs38DH-extra.fa > hs38DH.fa
            rm GCA_000001405.15_GRCh38_full_analysis_set.fna
        else
            {input}/run-gen-ref {wildcards.ref}
        fi
        
        if [ {params.region} != "all" ]; then
            samtools faidx {wildcards.ref}.fa {params.region} > {output}
            rm -f {wildcards.ref}.fa*
        else
            mv {wildcards.ref}.fa {output}
        fi
        """


# TODO: edit fix_names.py to also change the names for hg38
rule fix_names_clean:
    input:
        rules.gen_ref.output,
    output:
        fa="resources/{ref}/genome.fa",
        fai="resources/{ref}/genome.fa.fai",
        chromsizes="resources/{ref}/genome.genome",
    log:
        "resources/{ref}/fix_names.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/fix_names.py"


rule get_eul1db:
    input:
        "resources/eul1db_SRIP.txt",
    output:
        "resources/eul1db/insertions.bed",
    conda:
        "../envs/env.yml"
    log:
        "resources/eul1db/get_eul1db.log",
    script:
        "../scripts/get_eul1db.py"


rule liftover:
    input:
        srip="resources/{db}/insertions.bed",
        fa=expand(rules.fix_names_clean.output.fa, ref=config["genome"]["build"]),
    output:
        "resources/{db}/insertions_lifted.bed",
    log:
        "resources/{db}/liftover.log",
    conda:
        "../envs/env.yml"
    params:
        ref=config["genome"]["build"],
        chain=config["chain"],
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        workflow/scripts/liftover.sh \
            {input.srip} {input.fa} {output} {wildcards.db} {params.ref} {params.chain}
        """


rule get_windows:
    input:
        rules.liftover.output,
    output:
        "resources/{db}/windows.csv",
    log:
        "resources/{db}/get_windows.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_windows.py"


rule get_rmsk:
    input:
        "resources/{ref}/genome.genome",
    output:
        rmsk="resources/{ref}/rmsk.out",
        ref_l1="resources/{ref}/reference_l1.csv",
    log:
        "resources/{ref}/rmsk.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_rmsk.py"
