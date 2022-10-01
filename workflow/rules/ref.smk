rule gen_ref:
    output:
        f"resources/{{ref}}/{gen_ref_basename}.fa",
    log:
        "resources/{ref}/gen_ref.log",
    conda:
        "../envs/env.yml"
    params:
        region=config["ref"]["region"],
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
            samtools faidx $TMP/{wildcards.ref}.fa {params.region} > {output}
        else
            mv $TMP/{wildcards.ref}.fa {output}
        fi
        """


# TODO: edit fix_names.py to also change the names for hg38
rule fix_names_clean:
    input:
        fa=rules.gen_ref.output,
    output:
        fa=f"resources/{{ref}}/{gen_ref_basename}_hg19names.fa",
        fai=f"resources/{{ref}}/{gen_ref_basename}_hg19names.fa.fai",
        chromsizes=f"resources/{{ref}}/{gen_ref_basename}_hg19names.genome",
    log:
        "resources/{ref}/fix_names.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/fix_names.py"


rule get_eul1db:
    input:
        genome=expand(
            f"resources/{{ref}}/{ref_basename}.genome", ref=config["ref"]["build"]
        ),
        eul1db="resources/eul1db_SRIP.txt",
    output:
        "resources/eul1db/windows.csv",
    conda:
        "../envs/env.yml"
    log:
        "resources/eul1db/get_eul1db.log",
    script:
        "../scripts/get_eul1db.py"


rule get_rmsk:
    input:
        f"resources/{{ref}}/{ref_basename}.genome",
    output:
        rmsk="resources/{ref}/rmsk.txt.gz",
        ref_l1="resources/{ref}/reference_l1.csv",
    log:
        "resources/{ref}/rmsk.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_rmsk.py"
