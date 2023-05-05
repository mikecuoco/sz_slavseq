rule get_features:
    input:
        bam=rules.sambamba_sort.output[0],
        bai=rules.sambamba_index.output[0],
    params:
        **config["get_features"],
    output:
        "{outdir}/results/model/get_features/{donor}/{sample}.pqt",
    log:
        "{outdir}/results/model/get_features/{donor}/{sample}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_features.py"


def get_labels_input(wildcards):
    sample_ids = samples.loc[wildcards.donor]["sample_id"].values
    features = expand(
        rules.get_features.output,
        sample=sample_ids,
        allow_missing=True,
    )
    knrgl = rules.get_donor_knrgl.output[0]
    return {"features": features, "knrgl": knrgl}


rule get_labels:
    input:
        unpack(get_labels_input),
        rmsk=rules.run_rmsk.output[0],
        sv_blacklist=HTTP.remote(
            "https://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed",
            keep_local=True,
            static=True,
        ),
        segdups=HTTP.remote(
            "https://cf.10xgenomics.com/supp/genome/GRCh38/segdups.bedpe",
            keep_local=True,
            static=True,
        ),
    params:
        **config["get_features"],
    output:
        "{outdir}/results/model/get_labels/{donor}.pqt",
    log:
        "{outdir}/results/model/get_labels/{donor}.log",
    conda:
        "../envs/features.yml"
    threads: 8
    script:
        "../scripts/get_labels.py"


rule fit:
    input:
        expand(
            rules.get_labels.output,
            donor="CommonBrain",
            allow_missing=True,
        ),
    output:
        model="{outdir}/results/model/train/model.pkl",
        features="{outdir}/results/model/train/features.txt",
    log:
        "{outdir}/results/model/train/train.log",
    threads: 24
    benchmark:
        "{outdir}/results/model/train/train.benchmark"
    conda:
        "../envs/model.yml"
    script:
        "../scripts/fit.py"


rule predict:
    input:
        model=rules.fit.output.model,
        data=rules.get_features.output,
        sv_blacklist=HTTP.remote(
            "https://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed",
            keep_local=True,
            static=True,
        ),
        segdups=HTTP.remote(
            "https://cf.10xgenomics.com/supp/genome/GRCh38/segdups.bedpe",
            keep_local=True,
            static=True,
        ),
    output:
        "{outdir}/results/model/predict/{donor}/{sample}.pqt",
    conda:
        "../envs/model.yml"
    log:
        "{outdir}/results/model/predict/{donor}/{sample}.log",
    script:
        "../scripts/predict.py"
