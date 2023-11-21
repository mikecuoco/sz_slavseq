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
        "{outdir}/results/{genome}/{region}/{donor}/{sample}.pqt",
    log:
        "{outdir}/results/{genome}/{region}/{donor}/{sample}.log",
    conda:
        "../envs/features.yml"
    params:
        lambda wc: config[f"{wc.region}_params"],
    script:
        "../scripts/make_regions.py"


rule local_background:
    input:
        pqt=rules.make_regions.output[0],
        bam=rules.sambamba_sort.output[0],
        bai=rules.sambamba_index.output[0],
    output:
        "{outdir}/results/{genome}/{region}/{donor}/{sample}_bg.pqt",
    log:
        "{outdir}/results/{genome}/{region}/{donor}/{sample}_bg.log",
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

    # TODO: replace with local_background, skipping for now to finish prelim analysis
    # current error in local_background:
    #     Traceback (most recent call last):
    #   File "/iblm/logglun02/mcuoco/workflows/sz_slavseq/.snakemake/scripts/tmp2kquc_3d.local_background.py", line 49, in <module>
    #     f = sw.features(w)
    #   File "/iblm/logglun02/mcuoco/workflows/sz_slavseq/workflow/rules/../scripts/pyslavseq/sliding_window.py", line 351, in features
    #     f["3end_gini"] = gini(np.array(l["3end"], dtype=np.float64))
    #   File "/iblm/logglun02/mcuoco/workflows/sz_slavseq/workflow/rules/../scripts/pyslavseq/sliding_window.py", line 71, in gini
    #     if np.amin(array) < 0:
    #   File "/iblm/logglun02/mcuoco/workflows/sz_slavseq/.snakemake/conda/65cff2bb81647bc9d297d6550e298d36_/lib/python3.10/site-packages/numpy/core/fromnumeric.py", line 2970, in amin
    #     return _wrapreduction(a, np.minimum, 'min', axis, None, out,
    #   File "/iblm/logglun02/mcuoco/workflows/sz_slavseq/.snakemake/conda/65cff2bb81647bc9d297d6550e298d36_/lib/python3.10/site-packages/numpy/core/fromnumeric.py", line 88, in _wrapreduction
    #     return ufunc.reduce(obj, axis, dtype, out, **passkwargs)
    # ValueError: zero-size array to reduction operation minimum which has no identity
    return expand(
        rules.make_regions.output[0],
        # rules.local_background.output,
        region=wildcards.region,
        sample=cells,
        allow_missing=True,
    )


rule join_donor_regions:
    input:
        cells=get_donor_regions,
        en_motif_pos=rules.en_motif.output.pos,
        en_motif_neg=rules.en_motif.output.neg,
    output:
        "{outdir}/results/{genome}/{region}/{donor}.pqt",
    log:
        "{outdir}/results/{genome}/{region}/{donor}.log",
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
        "{outdir}/results/{genome}/{region}/{donor}_labelled.pqt",
    log:
        "{outdir}/results/{genome}/{region}/{donor}_labelled.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/label_donor_regions.py"


# TODO: add rule to generate report on regions
# By cells, donors, and labels
# 1. # regions
# 2. if peaks, length
# 3. reads / region
