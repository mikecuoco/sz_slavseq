rule get_hg38:
    output:
        expand("resources/GRCh38_full_analysis_set_plus_decoy_hla{ext}", ext=extensions),
        "resources/hg38.genome"
    log:
        "results/hg38.log"
    conda:
        "../envs/hg38.yml"
    shell:
        "workflow/scripts/hg38.sh > {log} 2>&1"