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
    "".join(config["genome"]["region"])
    if isinstance(config["genome"]["region"], list)
    else config["genome"]["region"]
)
region_name = f"_{region}" if region != "all" else ""
genome_name = config["genome"]["name"] + region_name


# generate hg38 reference with decoy and alt contigs
rule gen_ref:
    input:
        fa=FTP.remote(
            config["genome"]["ftp"],
            keep_local=True,
            static=True,
            immediate_close=True,
        ),
    output:
        multiext(
            f"{{outdir}}/resources/{genome_name}",
            ".fa",
            ".fa.fai",
            ".genome",
        ),
    log:
        "{outdir}/resources/gen_ref.log",
    conda:
        "../envs/ref.yml"
    params:
        region=" ".join(config["genome"]["region"])
        if isinstance(config["genome"]["region"], list)
        else config["genome"]["region"],
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
        multiext(
            f"{{outdir}}/resources/{genome_name}.fa",
            ".out",
            ".masked",
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
        speed="-s" if config["genome"]["region"] == "all" else "-qq",
    shadow:
        "shallow"
    threads: 24
    shell:
        """
        RepeatMasker -pa {threads} {params.speed} -lib {input.lib} -no_is -dir $(dirname {input.fa}) {input.fa} > {log} 2>&1
        """


rule get_donor_knrgl:
    input:
        lambda wc: donors.loc[wc.donor]["KNRGL"],
    output:
        "{outdir}/resources/{donor}_insertions.bed",
    conda:
        "../envs/ref.yml"
    log:
        "{outdir}/resources/{donor}/get_vcf.log",
    shell:
        """
        bcftools query -f "%CHROM\t%POS\t%END\t%STRAND\n" {input} | awk -v OFS='\t' '{{print $1,$2-1,$3,$4}}' > {output}
        """
