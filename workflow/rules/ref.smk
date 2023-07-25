# handle specified region
region = (
    "".join(config["genome"]["region"])
    if isinstance(config["genome"]["region"], list)
    else config["genome"]["region"]
)
region_name = f"_{region}" if region != "all" else ""
genome_name = config["genome"]["name"] + region_name
assert (
    "38" in genome_name
), "Only GRCh38 is supported due to segdup and blacklist regions used from 10x genomics"


rule get_genome:
    input:
        fa=FTP.remote(
            config["genome"]["fasta"],
            keep_local=True,
            static=True,
            immediate_close=True,
        ),
    output:
        fa=f"{{outdir}}/resources/{genome_name}.fa",
        fai=f"{{outdir}}/resources/{genome_name}.fa.fai",
    log:
        "{outdir}/resources/gen_ref.log",
    conda:
        "../envs/ref.yml"
    script:
        "../scripts/get_genome.py"


rule make_dfam_lib:
    output:
        "{outdir}/resources/LINE1_lib.fa",
    log:
        "{outdir}/resources/rmsk_lib.log",
    conda:
        "../envs/ref.yml"
    params:
        accessions=[
            "DF0000225",
            "DF0000339",
            "DF0000340",
            "DF0000341",
            "DF0000342",
            "DF0000343",
        ],
    shell:
        """
        touch {log} && exec 1>{log} 2>&1
        touch {output}
        for a in {params.accessions}; do
            curl -s https://dfam.org/api/families/$a | jq -r '.name' | awk '{{print ">"$1}}'  >> {output}
            curl -s https://dfam.org/api/families/$a | jq -r '.consensus_sequence' >> {output}
        done
        """


rule rmsk_to_bed:
    input:
        FTP.remote(
            config["genome"]["rmsk"],
            keep_local=True,
            static=True,
            immediate_close=True,
        ),
    output:
        rmsk="{outdir}/resources/rmsk.bed",
        rmsk_1kb_3end="{outdir}/resources/rmsk_1kb_3end.bed",
        rmsk_20kb="{outdir}/resources/rmsk_20kb.bed",
    log:
        "{outdir}/resources/rmsk_to_bed.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/rmsk_to_bed.py"


# get vcf, convert to bed, remove orphan insertions (higher FP rate and won't be detected in SLAV-seq)
rule xtea_to_bed:
    input:
        lambda wc: donors.loc[wc.donor]["xtea"],
    output:
        xtea="{outdir}/resources/{donor}_insertions.bed",
        xtea_1kb_3end="{outdir}/resources/{donor}_insertions_1kb_3end.bed",
        xtea_20kb="{outdir}/resources/{donor}_insertions_20kb.bed",
    conda:
        "../envs/features.yml"
    log:
        "{outdir}/resources/{donor}_to_bed.log",
    script:
        "../scripts/xtea_to_bed.py"
