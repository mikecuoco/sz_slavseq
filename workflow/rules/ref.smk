rule get_ref:
    output:
        multiext("resources/{ref}", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".genome")
    log:
        "resources/{ref}.log"
    conda:
        "../envs/env.yml"
    script:
        "../scripts/{wildcards.ref}.sh"