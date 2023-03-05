rule install_bwakit:
    output:
        directory("resources/bwa.kit"),
    conda:
        "../envs/ref.yml"
    log:
        "resources/install_bwakit.log",
    shell:
        """
        mkdir -p $(dirname {output}) && cd $(dirname {output})
        wget -O- -q --no-config https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2 | tar xfj -
        """


# handle specified region
region = (
    "".join(config["region"])
    if isinstance(config["region"], list)
    else config["region"]
)
region_name = f"_{region}" if region != "all" else ""

from snakemake.remote import FTP

FTP = FTP.RemoteProvider()


# generate hg38 reference with decoy and alt contigs
rule gen_ref:
    input:
        fa=FTP.remote(
            "ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
            keep_local=True,
            static=True,
        ),
    output:
        multiext(
            f"{{outdir}}/resources/hs38DH{region_name}",
            ".fa",
            ".fa.fai",
            ".genome",
        ),
    log:
        "{outdir}/resources/gen_ref.log",
    conda:
        "../envs/ref.yml"
    params:
        region=" ".join(config["region"])
        if isinstance(config["region"], list)
        else config["region"],
    shadow:
        "shallow"
    shell:
        """
        # start logging
        touch {log} && exec 2>{log}

        # filter for the region if specified
        if [ "{params.region}" != "all" ]; then
            samtools faidx {input} {params.region} > {output[0]}
            rm -f {input}
        else
            mv {input} {output[0]}
        fi

        # index
        samtools faidx {output[0]}
        cut -f 1,2 {output[1]} > {output[2]}
        """


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


rule run_rmsk:
    input:
        fa=rules.gen_ref.output[0],
        lib=rules.make_dfam_lib.output,
    output:
        protected(
            multiext(
                f"{{outdir}}/resources/hs38DH{region_name}.fa",
                ".out",
                ".masked",
            )
        ),
    log:
        "{outdir}/resources/run_rmsk.log",
    conda:
        "../envs/ref.yml"
    params:
        # -s Slow search; 0-5% more sensitive, 2-3 times slower than default;
        # empty string is default
        # -q Quick search; 5-10% less sensitive, 2-5 times faster than default
        # -qq Rush job; about 10% less sensitive, 4->10 times faster than default
        speed="-s" if config["region"] == "all" else "-qq",
    shadow:
        "shallow"
    threads: 24
    shell:
        """
        RepeatMasker -pa {threads} {params.speed} -lib {input.lib} -no_is -dir $(dirname {input.fa}) {input.fa} > {log} 2>&1
        """


rule get_eul1db:
    input:
        "resources/eul1db_SRIP.txt",
    output:
        "{outdir}/resources/eul1db/hg19_insertions.bed",
    conda:
        "../envs/ref.yml"
    log:
        "{outdir}/resources/get_eul1db.log",
    script:
        "../scripts/get_eul1db.py"


rule get_dbvar:
    output:
        vcf="{outdir}/resources/dbVar/GRCh38.variant_call.all.vcf.gz",
        tbi="{outdir}/resources/dbVar/GRCh38.variant_call.all.vcf.gz.tbi",
        bed="{outdir}/resources/dbVar/hs38DH_insertions.bed",
    conda:
        "../envs/ref.yml"
    log:
        "{outdir}/resources/get_dbvar.log",
    shell:
        """
        touch {log} && exec 1>{log} 2>&1
        curl -s https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.all.vcf.gz > {output.vcf}
        curl -s https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.all.vcf.gz.tbi > {output.tbi}

        bcftools query -f "%CHROM\t%POS\t%END\t%ALT\n" {output.vcf} | \
            grep "INS:ME:LINE1" | \
            uniq -u | \
            sed -e 's/^/chr/' | \
            awk -v OFS='\t' '{{print $1,$2,$3}}' > {output.bed}
        """


def get_KNRGL_build(wildcards):
    return config["KNRGL"][wildcards.db]["build"]


def get_liftover_input(wildcards):
    KNRGL_build = get_KNRGL_build(wildcards)
    return f"{wildcards.outdir}/resources/{wildcards.db}/{KNRGL_build}_insertions.bed"


rule liftover:
    input:
        get_liftover_input,
    output:
        "{outdir}/resources/{db}/{target}_lifted_insertions.bed",
    log:
        "{outdir}/resources/{db}/{target}_liftover.log",
    conda:
        "../envs/ref.yml"
    params:
        source=get_KNRGL_build,
    script:
        "../scripts/liftover_bed.sh"
