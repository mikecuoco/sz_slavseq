rule get_hg38:
    output:
        expand("resources/GRCh38_full_analysis_set_plus_decoy_hla.{ext}", ext=config["extensions])
    shell:
        "scripts/hg38.sh"