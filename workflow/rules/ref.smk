rule get_line1_consensus:
    output:
        hmm="resources/LINE1_3end/LINE1_3end.hmm",
        fa="resources/LINE1_3end/LINE1_3end.fa",
    log:
        "resources/LINE1_3end/dfam_query.log",
    conda:
        "../envs/ref.lock.yml"
    params:
        accessions=[
            "DF000000225",
            "DF000000339",
            "DF000000340",
            "DF000000341",
            "DF000000342",
            "DF000000343",
            "DF000000344",
            "DF000000345",
            "DF000000347",
            "DF000000327",
            "DF000000329",
            "DF000000331",
            "DF000000333",
            "DF000000335",
            "DF000000336",
            "DF000000337",
        ],
    shell:
        """
        touch {output} {log}
        accessions=({params.accessions})
        for a in "${{accessions[@]}}"; do
            curl -s https://dfam.org/api/families/$a/hmm?format=hmm >> {output.hmm} 2>> {log}
            curl -s https://dfam.org/api/families/$a | jq -r '.name' | awk '{{print ">"$1}}' >> {output.fa} 2>> {log}
            curl -s https://dfam.org/api/families/$a | jq -r '.consensus_sequence' >> {output.fa} 2>> {log}
        done
        """


# generate high accuracy annotation of L1 sequences in references
rule run_rmsk:
    input:
        fa=config["genome"]["fasta"],
        lib=rules.get_line1_consensus.output.hmm,
    output:
        out=multiext(
            config["genome"]["fasta"],
            ".out",
            ".masked",
        ),
        bed=config["genome"]["fasta"] + ".rmsk.bed",
    log:
        config["genome"]["fasta"] + ".rmsk.log",
    conda:
        "../envs/ref.lock.yml"
    params:
        # -s Slow search; 0-5% more sensitive, 2-3 times slower than default;
        # empty string is default
        # -q Quick search; 5-10% less sensitive, 2-5 times faster than default
        # -qq Rush job; about 10% less sensitive, 4->10 times faster than default
        speed="-qq" if config["genome"]["name"] == "chr2122" else "-s",
    shadow:
        "shallow"
    threads: 24
    shell:
        """
        RepeatMasker -pa {threads} -lib {input.lib} -no_is -e hmmer {params.speed} {input} > {log} 2>&1
        rmsk2bed < {output.out[0]} | cut -f 1-15 > {output.bed}
        """


rmsk_anno = {
    "l1hs": config["genome"]["fasta"] + ".rmsk.l1hs.bed",
    "l1pa2": config["genome"]["fasta"] + ".rmsk.l1pa2.bed",
    "l1pa3": config["genome"]["fasta"] + ".rmsk.l1pa3.bed",
    "l1pa4": config["genome"]["fasta"] + ".rmsk.l1pa4.bed",
    "l1pa5": config["genome"]["fasta"] + ".rmsk.l1pa5.bed",
    "l1pa6": config["genome"]["fasta"] + ".rmsk.l1pa6.bed",
    "polyA": config["genome"]["fasta"] + ".rmsk.polyA.bed",
    "polyT": config["genome"]["fasta"] + ".rmsk.polyT.bed",
}


rule filter_rmsk:
    input:
        rules.run_rmsk.output.bed,
    output:
        **rmsk_anno,
    log:
        config["genome"]["fasta"] + ".rmsk.filter.log",
    conda:
        "../envs/model.yml"
    script:
        "../scripts/filter_rmsk.py"


rule blast_primers:
    input:
        ref_fa=config["genome"]["fasta"],
        primer_fa="resources/L1_capture/L1_capture.fa",
    output:
        db=multiext(
            "resources/{genome}/BLASTDB/{genome}",
            ".ndb",
            ".nhr",
            ".nin",
            ".njs",
            ".not",
            ".nsq",
            ".ntf",
            ".nto",
        ),
        bed="resources/{genome}/blast_primers/L1_primers.bed",
    log:
        "resources/{genome}/blast_primers/blast_primers.log",
    conda:
        "../envs/ref.lock.yml"
    script:
        "../scripts/blast_primers.py"


rule en_motif:
    input:
        config["genome"]["fasta"],
    output:
        pos="resources/{genome}/en_pos_score.wig",
        neg="resources/{genome}/en_neg_score.wig",
    log:
        "resources/{genome}/en_motif.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/en_motif.py"


def get_vcf(wildcards):
    if "breakpoints" in wildcards.vcf:
        return donors.loc[wildcards.donor, "breakpoints"]
    else:
        return config["genome"][wildcards.vcf]


rule vcf2bed:
    input:
        vcf=get_vcf,
        meta=config["donors"],
        chrom_sizes=config["genome"]["genome"],
    output:
        "{outdir}/results/{genome}/vcf2bed/{donor}/{vcf}.bed",
    log:
        "{outdir}/results/{genome}/vcf2bed/{donor}/{vcf}.log",
    conda:
        "../envs/ref.lock.yml"
    params:
        slop_size=100,
    shell:
        """
        exec &> {log}

        LIBDID=$(csvcut -t -c "donor_id","libd_id" {input.meta} | csvgrep -c 1 -r "^{wildcards.donor}$" | csvcut -c 2 | tail -n +2)
        # check if the vcf is a xtea vcf
        if [ "{wildcards.vcf}" == "xtea" ]; then
            bcftools view -s "$LIBDID.md" {input.vcf} | \
                bcftools view -i "AC>0" | \
                bcftools query -f '%CHROM\t%POS\t%END\t%FILTER\t%SUBTYPE\t%AC\n' -e 'INFO/SUBTYPE ~ "transduction"' | \
                awk -v OFS='\t' '{{print $1, $2-1, $3, $4, $5, $6}}' | \
                bedtools slop -i - -g {input.chrom_sizes} -b {params.slop_size} > {output}
        elif [ "{wildcards.vcf}" == "megane_percentile" ] || [ "{wildcards.vcf}" == "megane_gaussian" ]; then
            bcftools view -s "$LIBDID" -i 'INFO/SVTYPE="LINE/L1"' {input.vcf} | \
                bcftools view -i "AC>0" | \
                bcftools query -f '%CHROM\t%0START\t%0END\t%FILTER\t%MEI\t%AC\n' | \
                bedtools slop -i - -g {input.chrom_sizes} -b {params.slop_size} > {output}
        elif [ "{wildcards.vcf}" == "megane_breakpoints" ]; then
            zcat {input.vcf} | \
                grep "L1HS" | \
                awk -v OFS='\t' '{{split($4,a,":"); split($5,b,":"); print $1, $2-1, $3, a[2], b[2], $6, $7}}' | \
                bedtools slop -i - -g {input.chrom_sizes} -b {params.slop_size} > {output}
        elif [ "{wildcards.vcf}" == "graffite" ]; then
            bcftools view -s "$LIBDID" -i 'INFO/repeat_ids="L1HS"' {input.vcf} | \
                bcftools view -i "AC>0" | \
                bcftools view -i 'strlen(ALT)>strlen(REF)' | \
                bcftools query -f '%CHROM\t%POS\t%END\t%FILTER\t%n_hits\t%total_match_span\n' | \
                bedtools slop -i - -g {input.chrom_sizes} -b {params.slop_size} > {output}
        fi
        """


def get_bed(wildcards):
    if "megane" in wildcards.vcf:
        return [
            *donors["megane_30x"].values,
            *donors["megane_90x_DLPFC"].dropna().values,
            *donors["megane_90x_HIPPO"].dropna().values,
        ]
    else:
        return expand(
            rules.vcf2bed.output, donor=donors["donor_id"], allow_missing=True
        )


rule merge_bed:
    input:
        get_bed,
    output:
        "{outdir}/results/{genome}/vcf2bed/{vcf}.merged.bed",
    log:
        "{outdir}/results/{genome}/vcf2bed/{vcf}.merge.log",
    conda:
        "../envs/ref.lock.yml"
    shell:
        """
        # if wildcards == megane
        if [[ {wildcards.vcf} == "megane" ]]; then
            cat {input} | grep "confidence:high;PASS" | grep "L1HS" | sort -k1,1 -k2,2n | bedtools merge -i - > {output} 2> {log}
        else
            cat {input} | sort -k1,1 -k2,2n | bedtools merge -i - > {output} 2> {log}
        fi
        """


knrgl_anno = {
    # "xtea": expand(
    # 	rules.merge_bed.output,
    # 	vcf="xtea",
    # 	allow_missing=True,
    # ),
    "megane": expand(
        rules.merge_bed.output,
        vcf="megane",
        allow_missing=True,
    ),
    "graffite": expand(
        rules.merge_bed.output,
        vcf="graffite",
        allow_missing=True,
    ),
}


rule annotation_report:
    input:
        **rmsk_anno,
        xtea=config["genome"]["xtea"],
        graffite=config["genome"]["graffite"],
        primer_sites=rules.blast_primers.output.bed,
    output:
        "{outdir}/results/{genome}/wgs_calls_report.ipynb",
    log:
        notebook="{outdir}/results/{genome}/wgs_calls_report.ipynb",
    conda:
        "../envs/model.yml"
    notebook:
        "../scripts/annotation_report.py.ipynb"


rule annotations:
    input:
        expand(
            rules.annotation_report.input,
            genome=config["genome"]["name"],
            outdir=config["outdir"],
        ),
