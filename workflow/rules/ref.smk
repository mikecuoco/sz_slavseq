rule gen_ref:
    output:
        multiext(f"resources/{{ref}}/{gen_ref_basename}", ".fa", ".fa.fai", ".genome"),
    log:
        "resources/{ref}/gen_ref.log",
    conda:
        "../envs/ref.yml"
    params:
        region=region,
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


rule run_rmsk:
    input:
        rules.gen_ref.output[0],
    output:
        multiext(f"resources/{{ref}}/{gen_ref_basename}.fa", ".out", ".masked"),
    log:
        "resources/{ref}/run_rmsk.log",
    conda:
        "../envs/ref.yml"
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
        touch {log} && exec 1>{log} 2>&1

        # download dfam
        wget -O- https://www.dfam.org/releases/Dfam_3.6/families/Dfam-p1_curatedonly.h5.gz | \
            gzip -dc > $CONDA_PREFIX/share/RepeatMasker/Libraries/Dfam.h5 

        # run RepeatMasker
        RepeatMasker -pa {threads} {params.speed} -species human -dir $(dirname {input}) {input} > {log} 2>&1

        # TODO: convert to bed
        """


rule get_eul1db:
    input:
        "resources/eul1db_SRIP.txt",
    output:
        "resources/hg19/hg19_eul1db_insertions.bed",
    conda:
        "../envs/ref.yml"
    log:
        "resources/get_eul1db.log",
    script:
        "../scripts/get_eul1db.py"


rule get_dbvar:
    output:
        vcf="resources/hs38DH/dbVar.variant_call.all.vcf.gz",
        tbi="resources/hs38DH/dbVar.variant_call.all.vcf.gz.tbi",
        bed="resources/hs38DH/hs38DH_dbVar_insertions.bed",
    conda:
        "../envs/ref.yml"
    log:
        "resources/get_dbvar.log",
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        URL=https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.all.vcf.gz
        wget -O- $URL > {output.vcf}
        wget -O- $URL.tbi > {output.tbi}
        bcftools query -f "%CHROM\t%POS\t%END\t%ALT\n" {output.vcf} | \
            grep "INS:ME:LINE1" | \
            uniq -u | \
            sed -e 's/^/chr/' | \
            awk '{{print $1,$2,$3}}' > {output.bed} 
        """


target_build = (
    "hg19" if config["genome"]["build"] == "hs37d5" else config["genome"]["build"]
)


rule liftover:
    input:
        get_liftover_input,
    output:
        expand("resources/{target}/{target}_{{db}}_insertions.bed", target=target_build),
    log:
        "resources/{db}_liftover.log",
    conda:
        "../envs/ref.yml"
    params:
        source=config["non_ref_germline_l1"]["build"],
        target=target_build,
    script:
        "../scripts/liftover_bed.sh"


rule fix_names:
    input:
        bed=get_fixnames_input,
        chrom_map="resources/hs37d5_map.tsv",
    output:
        "resources/{ref}/{ref}_{db}_insertions.bed",
    log:
        "resources/{ref}/{ref}_{db}_fixnames.log",
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
