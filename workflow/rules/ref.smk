rule gen_ref:
    output:
        multiext(f"resources/{{ref}}/{gen_ref_basename}", ".fa", ".fa.fai", ".genome"),
    log:
        "resources/{ref}/gen_ref.log",
    conda:
        "../envs/env.yml"
    params:
        region=" ".join(config["genome"]["region"])
        if isinstance(config["genome"]["region"], list)
        else config["genome"]["region"],
    shadow:
        "shallow"
    cache: True
    shell:
        """
        # start logging
        touch {log} && exec 2>{log} 

        # allow errors temporarily to handle broken pipe with hs37d5
        set +e 
        bash workflow/scripts/run-gen-ref.sh {wildcards.ref}
        set -e

        # filter for the region if specified
        if [ "{params.region}" != "all" ]; then
            samtools faidx {wildcards.ref}.fa {params.region} > {output[0]}
        else
            mv {wildcards.ref}.fa {output[0]}
        fi

        samtools faidx {output[0]}
        cut -f 1,2 {output[1]} > {output[2]} # TODO: use samtools dict instead?
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
        multiext("resources/{ref}/{ref}_{db}_insertions", ".bed", ".bed.unmap"),  # is unmap file always generated?
    log:
        "resources/{ref}/{ref}_{db}_liftover.log",
    conda:
        "../envs/env.yml"
    params:
        source_build=config["non_ref_germline_l1"]["build"],
        target_build=config["genome"]["build"],  # this is the same as wildcards.ref
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        bash workflow/scripts/liftover_bed.sh -s {params.source_build} -t {params.target_build} -i {input} -o {output[0]} 
        """


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
        speed="-s" if config["genome"]["region"] == "all" else "-qq",
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
