rule install_bwakit:
    output:
        directory("resources/bwa.kit"),
    conda:
        "../envs/env.yml"
    log:
        "resources/install_bwakit.log",
    shell:
        """
        mkdir -p resources && cd resources
        wget -O- -q --no-config https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2 | tar xfj -
        """


rule gen_ref:
    input:
        rules.install_bwakit.output,
    output:
        "resources/{ref}/genome_og.fa",
    log:
        "resources/{ref}/gen_ref.log",
    conda:
        "../envs/env.yml"
    params:
        region=config["genome"]["region"],
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        # run bwa.kit function
        {input}/run-gen-ref {wildcards.ref}
        if [ {params.region} != "all" ]; then
            samtools faidx {wildcards.ref}.fa {params.region} > {output}
            rm -f {wildcards.ref}.fa*
        else
            mv {wildcards.ref}.fa {output}
        fi
        """


# TODO: edit fix_names.py to also change the names for hg38
rule fix_names_clean:
    input:
        rules.gen_ref.output,
    output:
        fa="resources/{ref}/genome.fa",
        fai="resources/{ref}/genome.fa.fai",
        chromsizes="resources/{ref}/genome.genome",
    log:
        "resources/{ref}/fix_names.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/fix_names.py"


rule get_eul1db:
    input:
        rules.fix_names_clean.output.chromsizes
    output:
        "resources/eul1db/insertions.bed",
        # "resources/eul1db/windows.csv",
    conda:
        "../envs/env.yml"
    log:
        "resources/eul1db/get_eul1db.log",
    script:
        "../scripts/get_eul1db.py"


rule liftover:
    input:
        insert=expand("resources/{db}/insertions.bed", db=config["ref"]["database"]),
        fa=expand(rules.fix_names_clean.output.fa, ref=config["genome"]["build"]),
    output:
        expand("resources/{db}/insertions_stable.vcf", db=config["ref"]["database"]),
    log:
        expand("resources/{db}/liftover.log", db=config["ref"]["database"]),
    conda:
        "../envs/env.yml"
    params:
        build=config["genome"]["build"]
        source=config["germline_line1"]["source"]
        # chain=config["chain"]
    shell:
        """
        touch {log} && exec 1>{log} 2>&1

        if [[ {params.chain} != None ]]; then
            wget -P resources/liftover/ -q --no-config \
                https://raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh37.novel_CUPs.bed

            Rscript workflow/scripts/bed_to_vcf.R {input.insert} resources/liftover/insertions.vcf \
                resources/{wildcards.ref}/genome.fa

            vcftools --vcf insertions.vcf --exclude-bed FASTA_BED.ALL_GRCh3N.novel_CUPs.bed \
                --recode --recode-INFO-all --out insertions_stable.vcf

            picard CreateSequenceDictionary -R {input.fa} -O resources/{wildcards.ref}/genome.dict

            picard LiftoverVcf -I insertions_stable.vcf -O insertions_lifted.vcf \
                -C resources/{params.chain} \
                --REJECT rejected_insertions.vcf -R {input.fa}
        else
            # make output consistent for all cases
            Rscript workflow/scripts/bed_to_vcf.R {input.insert} resources/liftover/insertions_stable.vcf \
                resources/{wildcards.ref}/genome.fa
        fi
        """
    

rule get_windows:
    input: 
        rules.liftover.output
    output:
        expand("resources/{db}/windows.csv", db=config["ref"]["database"]),
    log:
        expand("resources/{db}/get_windows.log", db=config["ref"]["database"]),
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_windows.py"


rule get_rmsk:
    input:
        "resources/{ref}/genome.genome",
    output:
        rmsk="resources/{ref}/rmsk.txt.gz",
        ref_l1="resources/{ref}/reference_l1.csv",
    log:
        "resources/{ref}/rmsk.log",
    conda:
        "../envs/env.yml"
    script:
        "../scripts/get_rmsk.py"
