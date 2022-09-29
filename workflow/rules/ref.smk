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
        "resources/{ref}/genome.fa",
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
            wget -O GRCh38_full_analysis_set.fna.gz -q --no-config \
                ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
            gunzip GRCh38_full_analysis_set.fna.gz
            cat GRCh38_full_analysis_set.fna {input}/resource-GRCh38/hs38DH-extra.fa > hs38DH.fa
            rm GRCh38_full_analysis_set.fna
        else
            wget -O hs37d5.fa.gz -q --no-config \
                ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
            gunzip hs37d5.fa.gz
        fi

        if [ {params.region} != "all" ]; then
            samtools faidx {wildcards.ref}.fa {params.region} > {output}
            rm -f {wildcards.ref}.fa*
        else
            mv {wildcards.ref}.fa {output}
        fi
        """


rule index_genome:
    input:
        rules.gen_ref.output,
    output:
        fai="resources/{ref}/genome.fa.fai",
        chromsizes="resources/{ref}/genome.genome",
    log:
        "resources/{ref}/index_genome.log",
    conda:
        "../envs/env.yml"
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        samtools faidx {input}
        cut -f 1,2 {output.fai} > {output.chromsizes}
        """


rule get_eul1db:
    input:
        "resources/eul1db_SRIP.txt",
    output:
        "resources/eul1db/insertions_hg19.bed",
    conda:
        "../envs/env.yml"
    log:
        "resources/eul1db/get_eul1db.log",
    script:
        "../scripts/get_eul1db.py"


rule fix_names:
    input:
        "resources/{db}/insertions_hg19.bed",
    output:
        "resources/{db}/insertions_hs37d5.bed",
    log:
        "resources/{db}/fix_names.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/fix_names.py"


rule liftover:
    input:
        non_ref_l1=get_non_ref_l1_for_liftover,
        fa=expand("resources/{ref}/genome.fa", ref=config["genome"]["build"]),
    output:
        "resources/{db}/insertions_lifted.bed",
    log:
        "resources/{db}/liftover.log",
    conda:
        "../envs/env.yml"
    params:
        chain=config["chain"],
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        if [[ {params.chain} == "hg19ToHg38.over.chain.gz" ]]; then
            # download chain file
            wget -O resources/hg19ToHg38.over.chain.gz -q --no-config \
                https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz 

            # download CUPs
            wget -O resources/GRCh37.novel_CUPs.bed -q --no-config \
                https://raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh37.novel_CUPs.bed

            grep -v -f resources/GRCh37.novel_CUPs.bed {input.non_ref_l1} > resources/{wildcards.db}/insertions_stable.bed
            CrossMap.py bed resources/{params.chain} resources/{wildcards.db}/insertions_stable.bed {output}
        else
            # no liftover necessary
            mv {input.non_ref_l1} {output}
        fi
        """


rule get_windows:
    input:
        srip=rules.liftover.output,
    output:
        "resources/{db}/windows.csv",
    log:
        "resources/{db}/get_windows.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_windows.py"


rule run_rmsk:
    input:
        "resources/{ref}/genome.fa",
    output:
        multiext("resources/{ref}/", "rmsk.out", "rmsk.align"),
    log:
        "resources/{ref}/run_rmsk.log",
    conda:
        "../envs/env.yml"
    cache: True
    threads: 16
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        # download dfam
        wget -O- https://www.dfam.org/releases/Dfam_3.6/families/Dfam-p1_curatedonly.h5.gz | \
            gzip -dc > $CONDA_PREFIX/share/RepeatMasker/Libraries/Dfam.h5 

        # run RepeatMasker
        RepeatMasker -pa {threads} -s -nolow -species human -dir $(dirname {input}) {input}
        """


rule get_rmsk:
    input:
        rmsk=rules.run_rmsk.output[0],
        genome="resources/{ref}/genome.genome",
    output:
        "resources/{ref}/reference_l1.csv",
    log:
        "resources/{ref}/get_rmsk.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_rmsk.py"
