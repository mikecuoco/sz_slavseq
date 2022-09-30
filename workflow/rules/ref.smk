rule gen_ref:
    output:
        "resources/{ref}_{region}/genome.fa",
    log:
        "resources/{ref}_{region}/gen_ref.log",
    conda:
        "../envs/env.yml"
    cache: True
    shell:
        """
        # run the script
        bash workflow/scripts/run-gen-ref.sh {wildcards.ref} 

        # filter for the region if specified
        if [ {wildcards.region} != "all" ]; then
            samtools faidx {wildcards.ref}.fa {wildcards.region} > {output}
            rm -f {wildcards.ref}.fa*
        else
            mv {wildcards.ref}.fa {output}
        fi
        """


# TODO: edit fix_names.py to also change the names for hg38
rule fix_names_clean:
    input:
        fa=rules.gen_ref.output,
    output:
        fa="resources/{ref}_{region}/genome.fa",
        fai="resources/{ref}_{region}/genome.fa.fai",
        chromsizes="resources/{ref}_{region}/genome.genome",
    log:
        "resources/{ref}_{region}/fix_names.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/fix_names.py"


rule get_eul1db:
    input:
        genome=expand(
            "resources/{ref}_{region}/genome.genome",
            ref=config["ref"]["build"],
            region=config["ref"]["region"],
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
        "resources/{ref}_{region}/genome.genome",
    output:
        rmsk="resources/{ref}_{region}/rmsk.txt.gz",
        ref_l1="resources/{ref}_{region}/reference_l1.csv",
    log:
        "resources/{ref}_{region}/rmsk.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_rmsk.py"
