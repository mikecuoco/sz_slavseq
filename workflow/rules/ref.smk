rule gen_ref:
    output:
        multiext(f"resources/{{ref}}/{gen_ref_basename}", ".fa", ".fa.fai", ".genome"),
    log:
        "resources/{ref}/gen_ref.log",
    conda:
        "../envs/env.yml"
    params:
        region=config["genome"]["region"],
    cache: True
    shell:
        """
        # start logging
        touch {log} && exec 2>{log} 

        # create temp dir, and clean up on exit
        TMP=$(mktemp -d)
        trap 'rm -rf -- "$TMP"' EXIT

        # run the script in the temp dir
        cd $TMP
        set +e # allow errors temporarily to handle broken pipe with hs37d5
        bash $OLDPWD/workflow/scripts/run-gen-ref.sh {wildcards.ref}
        set -e

        # filter for the region if specified
        cd $OLDPWD
        if [ {params.region} != "all" ]; then
            samtools faidx $TMP/{wildcards.ref}.fa {params.region} > {output[0]}
        else
            mv $TMP/{wildcards.ref}.fa {output[0]}
        fi

        samtools faidx {output[0]}
        cut -f 1,2 {output[1]} > {output[2]}
        """


rule get_eul1db:
    input:
        "resources/eul1db_SRIP.txt",
    output:
        non_ref_l1_bed,
    conda:
        "../envs/env.yml"
    params:
        ref=config["genome"]["build"],
    log:
        "resources/get_eul1db.log",
    script:
        "../scripts/get_eul1db.py"


rule liftover:
    input:
        non_ref_l1_bed,
    output:
        multiext(non_ref_l1_bed_final, "", ".unmap"),
    log:
        "resources/{ref}/{ref}_{db}_liftover.log",
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

            CrossMap.py bed resources/{params.chain} <(grep -v -f resources/GRCh37.novel_CUPs.bed {input}) {output[0]}
        else
            echo "unable to liftover, chain file unavailable" && exit 1
        fi
        """


rule get_non_ref_l1_windows:
    input:
        non_ref_l1=non_ref_l1_bed_final,
        chromsizes=rules.gen_ref.output[2],
    output:
        "resources/{ref}/{ref}_{db}_windows.csv",
    log:
        "resources/{ref}/{ref}_{db}_windows.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_windows.py"


rule run_rmsk:
    input:
        rules.gen_ref.output[0],
    output:
        multiext(f"resources/{{ref}}/{gen_ref_basename}.fa", ".out", ".masked"),
    log:
        "resources/{ref}/run_rmsk.log",
    conda:
        "../envs/env.yml"
    params:
        # -s Slow search; 0-5% more sensitive, 2-3 times slower than default;
        # empty string is default
        # -q Quick search; 5-10% less sensitive, 2-5 times faster than default
        # -qq Rush job; about 10% less sensitive, 4->10 times faster than default
        speed="-qq",
    cache: True
    threads: 16
    shell:
        """
        # download dfam
        wget -O- https://www.dfam.org/releases/Dfam_3.6/families/Dfam-p1_curatedonly.h5.gz | \
            gzip -dc > $CONDA_PREFIX/share/RepeatMasker/Libraries/Dfam.h5 

        # run RepeatMasker
        RepeatMasker -pa {threads} {params.speed} -species human -dir $(dirname {input}) {input} > {log} 2>&1

        # TODO: convert to bed
        """


rule get_rmsk_windows:
    input:
        rmsk=rules.run_rmsk.output[0],
        chromsizes=rules.gen_ref.output[2],
    output:
        "resources/{ref}/reference_l1.csv",
    log:
        "resources/{ref}/get_rmsk.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_rmsk.py"
