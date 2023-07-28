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
    }


rule call_bulk_peaks:
    input:
        unpack(get_bulk_sample),
    output:
        pqt="{outdir}/results/bulk_peaks/{donor}.pqt",
        bed="{outdir}/results/bulk_peaks/{donor}.bed",
    log:
        "{outdir}/results/bulk_peaks/{donor}.log",
    params:
        **config["get_features"],
    conda:
        "../envs/model.yml"
    script:
        "../scripts/call_bulk_peaks.py"


rule get_features:
    input:
        bam=rules.sambamba_sort.output[0],
        bai=rules.sambamba_index.output[0],
    output:
        windows="{outdir}/results/model/get_features/{donor}/{sample}_windows.pqt",
        # peaks="{outdir}/results/model/get_features/{donor}/{sample}_peaks.pqt",
    log:
        "{outdir}/results/model/get_features/{donor}/{sample}.log",
    params:
        **config["get_features"],
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_features.py"


def get_donor_features(wildcards):
    donor_samples = samples.loc[samples["donor_id"] == wildcards.donor]
    cells = donor_samples[~donor_samples["sample_id"].str.contains("gDNA")]["sample_id"]

    res = {
        "features": expand(
            rules.get_features.output.windows,
            sample=cells,
            allow_missing=True,
        ),
    }
    # add bulk peaks if not common brain
    if wildcards.donor != "CommonBrain":
        res["bulk_peaks"] = rules.call_bulk_peaks.output.pqt

    return res


rule get_labels:
    input:
        unpack(get_donor_features),
        xtea=rules.xtea_to_bed.output.xtea,
        xtea_1kb_3end=rules.xtea_to_bed.output.xtea_1kb_3end,
        xtea_20kb=rules.xtea_to_bed.output.xtea_20kb,
        rmsk=rules.rmsk_to_bed.output.rmsk,
        rmsk_1kb_3end=rules.rmsk_to_bed.output.rmsk_1kb_3end,
        rmsk_20kb=rules.rmsk_to_bed.output.rmsk_20kb,
    params:
        **config["get_features"],
    output:
        "{outdir}/results/model/get_labels/{donor}.pqt",
    log:
        "{outdir}/results/model/get_labels/{donor}.log",
    conda:
        "../envs/model.yml"
    threads: 16
    script:
        "../scripts/get_labels.py"


rule fit:
    input:
        expand(
            rules.get_labels.output,
            donor=donors.loc[~donors["xtea"].isnull(), "donor_id"],
            outdir=config["outdir"],
        ),
    output:
        model="{outdir}/results/model/train/model.pkl",
        features="{outdir}/results/model/train/features.txt",
    log:
        "{outdir}/results/model/train/train.log",
    threads: 24
    params:
        # bad÷_cells
        chromosomes=["chr{}".format(i) for i in range(1, 22)],
    conda:
        "../envs/model.yml"
    script:
        "../scripts/fit.py"
