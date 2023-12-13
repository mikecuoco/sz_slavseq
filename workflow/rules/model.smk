def get_bulk_sample(wildcards):
    donor_samples = samples.loc[samples["donor_id"] == wildcards.donor]
    bulk = donor_samples[donor_samples["sample_id"].str.contains("gDNA")]["sample_id"]
    return {
        "bam": expand(
            rules.sambamba_sort.output[0],
            sample=bulk,
            allow_missing=True,
        ),
        "bai": expand(
            rules.sambamba_index.output[0],
            sample=bulk,
            allow_missing=True,
        ),
        "primer_bed": rules.blast_primers.output.bed,
    }


rule call_bulk_peaks:
    input:
        unpack(get_bulk_sample),
    output:
        pqt="{outdir}/results/{genome}/bulk_peaks/{donor}.pqt",
        bed="{outdir}/results/{genome}/bulk_peaks/{donor}.bed",
    log:
        "{outdir}/results/{genome}/bulk_peaks/{donor}.log",
    conda:
        "../envs/model.yml"
    params:
        **config["bulk_peaks_params"],
    script:
        "../scripts/call_bulk_peaks.py"


rule make_regions:
    input:
        bam=rules.sambamba_sort.output[0],
        bai=rules.sambamba_index.output[0],
    output:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}.pqt",
    log:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}.log",
    conda:
        "../envs/features.yml"
    params:
        regions_params.instance,
    script:
        "../scripts/make_regions.py"


rule local_background:
    input:
        pqt=rules.make_regions.output[0],
        bam=rules.sambamba_sort.output[0],
        bai=rules.sambamba_index.output[0],
    output:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}_bg.pqt",
    log:
        "{outdir}/results/{genome}/{params}/{donor}/{sample}_bg.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/local_background.py"


with open("resources/bad_cells.txt", "r") as f:
    bad_cells = [line.strip() for line in f.readlines()]


def get_donor_regions(wildcards):
    donor_cells = samples.loc[
        (samples["donor_id"] == wildcards.donor)
        & (~samples["sample_id"].str.contains("gDNA"))
    ]["sample_id"].values
    cells = [c for c in donor_cells if c not in bad_cells]

    return expand(
        rules.make_regions.output[0],
        # rules.local_background.output,
        sample=cells,
        allow_missing=True,
    )


rule join_donor_regions:
    input:
        cells=get_donor_regions,
        en_motif_pos=rules.en_motif.output.pos,
        en_motif_neg=rules.en_motif.output.neg,
    output:
        "{outdir}/results/{genome}/{params}/{donor}.pqt",
    log:
        "{outdir}/results/{genome}/{params}/{donor}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/join_donor_regions.py"


# TODO: add genome-specific mappability regions from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3 (specificy in config)
rule label_donor_regions:
    input:
        regions=rules.join_donor_regions.output,
        bulk_peaks=rules.call_bulk_peaks.output.bed,
        primer_sites=rules.blast_primers.output.bed,
        xtea_vcf=config["genome"]["xtea"],
        megane_vcf=config["genome"]["megane"],
        rmsk=rules.run_rmsk.output[0],
    output:
        "{outdir}/results/{genome}/{params}/{donor}_labelled.pqt",
    log:
        "{outdir}/results/{genome}/{params}/{donor}_labelled.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/label_donor_regions.py"


# TODO: add rule to generate report on regions
# By cells, donors, and labels
# 1. # regions
# 2. if peaks, length
# 3. reads / region
