rule get_ref:
    output:
        multiext("resources/{ref}", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".genome")
    log:
        "resources/{ref}.log"
    conda:
        "../envs/env.yml"
    shell:
        "if [[ {wildcards.ref} == 'GRCh38' ]]; then workflow/scripts/GRCh38.sh > {log} 2>&1; else workflow/scripts/hs37d5.sh > {log} 2>&1; fi"