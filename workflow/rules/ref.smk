# handle specified region
region = (
    "".join(config["genome"]["region"])
    if isinstance(config["genome"]["region"], list)
    else config["genome"]["region"]
)
region_name = f"_{region}" if region != "all" else ""

from snakemake.remote import FTP

FTP = FTP.RemoteProvider()


rule gen_ref:
    output:
        multiext(
            f"{{outdir}}/resources/{{ref}}/{{ref}}{region_name}",
            ".fa",
            ".fa.fai",
            ".genome",
        ),
    log:
        "{outdir}/resources/{ref}/gen_ref.log",
    conda:
        "../envs/ref.yml"
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

        # hs37d5 only accepts strings of digits (e.g. '22')
        # delete letters from params.region if necessary (e.g. 'chr22' -> '22')
        if [ "{params.region}" != "all" ] && [ {wildcards.ref} == "hs37d5" ]; then
            region=$(echo "{params.region}" | sed 's/[a-z]//gI')
        else
            region="{params.region}"
        fi

        # filter for the region if specified
        if [ "{params.region}" != "all" ]; then
            samtools faidx {wildcards.ref}.fa $region > {output[0]}
        else
            mv {wildcards.ref}.fa {output[0]}
        fi

        samtools faidx {output[0]}
        cut -f 1,2 {output[1]} > {output[2]} # TODO: use samtools dict instead?
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
            f"{{outdir}}/resources/{{ref}}/{{ref}}{region_name}.fa",
            ".out",
            ".masked",
        ),
    log:
        "{outdir}/resources/{ref}/run_rmsk.log",
    conda:
        "../envs/ref.yml"
    params:
        # -s Slow search; 0-5% more sensitive, 2-3 times slower than default;
        # empty string is default
        # -q Quick search; 5-10% less sensitive, 2-5 times faster than default
        # -qq Rush job; about 10% less sensitive, 4->10 times faster than default
        speed="-s" if config["genome"]["region"] == "all" else "-qq",
    cache: True
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
    input:
        FTP.remote(
            [
                "ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.all.vcf.gz",
                "ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.all.vcf.gz.tbi",
            ],
            static=True,
        ),
    output:
        "{outdir}/resources/dbVar/hs38DH_insertions.bed",
    conda:
        "../envs/ref.yml"
    log:
        "{outdir}/resources/get_dbvar.log",
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        bcftools query -f "%CHROM\t%POS\t%END\t%ALT\n" {input[0]} | \
            grep "INS:ME:LINE1" | \
            uniq -u | \
            sed -e 's/^/chr/' | \
            awk -v OFS='\t' '{{print $1,$2,$3}}' > {output} 
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


def get_fixnames_input(wildcards):
    KNRGL_build = get_KNRGL_build(wildcards)
    if wildcards.ref == "hs37d5" and KNRGL_build == "hg19":
        return f"{wildcards.outdir}/resources/{wildcards.db}/hg19_insertions.bed"
    elif wildcards.ref == "hs37d5" and KNRGL_build != "hg19":
        return f"{wildcards.outdir}/resources/{wildcards.db}/hg19_lifted_insertions.bed"


rule fix_names:
    input:
        bed=get_fixnames_input,
        chrom_map="resources/hs37d5_map.tsv",
    output:
        "{outdir}/resources/{db}/{ref}_fixnames_insertions.bed",
    log:
        "{outdir}/resources/{db}/{ref}_fixnames.log",
    run:
        # read in the bed file
        bed = pd.read_csv(input["bed"], sep="\t", names=["chr", "start", "end"])
        # read in the chromosome map
        chrom_map = pd.read_csv(input["chrom_map"], sep="\t", names=["hs37d5", "ann"])
        # change the names in a loop
        for name in chrom_map["ann"].to_list():
            bed.loc[bed["chr"] == name, "chr"] = chrom_map.loc[
                chrom_map["ann"] == name, "hs37d5"
            ].values[0]

        bed.to_csv(output[0], sep="\t", index=False, header=False)
