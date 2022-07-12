rule get_ref:
    output:
        expand("resources/{ref}{ext}", ext=extensions, ref=config["ref"])
    log:
        expand("results/{ref}.log", ref=config["ref"])
    conda:
        "../envs/ref.yml"
    shell:
        "if [[ {config[ref]} == 'GRCh38' ]]; then workflow/scripts/GRCh38.sh > {log} 2>&1; else workflow/scripts/hs37d5.sh > {log} 2>&1; fi"