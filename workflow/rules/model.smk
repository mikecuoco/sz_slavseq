rule get_features:
    input:
        bam=rules.sort.output[0],
        bai=rules.index.output[0],
        fa=rules.gen_ref.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        **config["get_features"],
    output:
        "{outdir}/results/model/get_features/{ref}_{db}/{donor}/{dna_type}/{sample}.pickle.gz",
    log:
        "{outdir}/results/model/get_features/{ref}_{db}/{donor}/{dna_type}/{sample}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_features.py"


def get_non_ref_l1(wildcards):
    KNRGL_build = get_KNRGL_build(wildcards)
    if wildcards.ref == "hs37d5":
        return f"{wildcards.outdir}/resources/{wildcards.db}/{wildcards.ref}_fixnames_insertions.bed"
    elif wildcards.ref != KNRGL_build:
        return f"{wildcards.outdir}/resources/{wildcards.db}/{wildcards.ref}_lifted_insertions.bed"
    else:
        return f"{wildcards.outdir}/resources/{wildcards.db}/{wildcards.ref}_insertions.bed"


def get_labels_input(wildcards):
    donor_samples = samples.loc[samples["donor"] == wildcards.donor]
    return {
        "bam": expand(
            rules.sort.output,
            sample=donor_samples.loc[samples["dna_type"] == "bulk"]["sample"],
            dna_type="bulk",
            allow_missing=True,
        ),
        "bai": expand(
            rules.index.output,
            sample=donor_samples.loc[samples["dna_type"] == "bulk"]["sample"],
            dna_type="bulk",
            allow_missing=True,
        ),
        "features": expand(
            rules.get_features.output,
            sample=donor_samples.loc[samples["dna_type"] == "mda"]["sample"],
            dna_type="mda",
            allow_missing=True,
        ),
    }


rule get_labels:
    input:
        unpack(get_labels_input),
        non_ref_l1=get_non_ref_l1,
        ref_l1=rules.run_rmsk.output[0],
        chromsizes=rules.gen_ref.output[2],
    params:
        **config["get_features"],
    output:
        "{outdir}/results/model/get_labels/{ref}_{db}/{donor}.pickle.gz",
    log:
        "{outdir}/results/model/get_labels/{ref}_{db}/{donor}.log",
    conda:
        "../envs/features.yml"
    script:
        "../scripts/get_labels.py"


rule folds:
    input:
        samples=expand(
            "{outdir}/results/model/get_labels/{ref}_{db}/{donor}.pickle.gz",
            donor=set(samples["donor"]),
            allow_missing=True,
        ),
    params:
        num_folds=config["num_folds"],
        min_reads=config["get_features"]["min_reads"],
    output:
        features="{outdir}/results/model/folds/{ref}_{db}/features.pickle",
        labels="{outdir}/results/model/folds/{ref}_{db}/labels.pickle",
        label_encoder="{outdir}/results/model/folds/{ref}_{db}/label_encoder.pickle",
        classes_per_cell="{outdir}/results/model/folds/{ref}_{db}/classes_per_cell.png",
        classes_per_donor="{outdir}/results/model/folds/{ref}_{db}/classes_per_donor.png",
        features_per_class="{outdir}/results/model/folds/{ref}_{db}/features_per_class.png",
        classes_per_donor_per_fold="{outdir}/results/model/folds/{ref}_{db}/classes_per_donor_per_fold.png",
    log:
        "{outdir}/results/model/folds/{ref}_{db}/folds.log",
    conda:
        "../envs/model.yml"
    wildcard_constraints:
        dna_type="\w+",
    script:
        "../scripts/folds.py"


rule train_test:
    input:
        features=rules.folds.output.features,
        labels=rules.folds.output.labels,
    params:
        model_params=lambda wc: config["models"][wc.model_id]["params"],
        model_name=lambda wc: config["models"][wc.model_id]["name"],
    output:
        # model="results/train_test/{ref}/{dna_type}/{model}/model.pickle",
        pred="{outdir}/results/model/train_test/{ref}_{db}/{model_id}/pred.pickle",
        proba="{outdir}/results/model/train_test/{ref}_{db}/{model_id}/proba.pickle",
    log:
        "{outdir}/results/model/train_test/{ref}_{db}/{model_id}.log",
    threads: 8
    benchmark:
        "{outdir}/results/model/train_test/{ref}_{db}/{model_id}.benchmark.txt"
    conda:
        "../envs/model.yml"
    script:
        "../scripts/train_test.py"


rule prcurve:
    input:
        labels=rules.folds.output.labels,
        proba=rules.train_test.output.proba,
        label_encoder=rules.folds.output.label_encoder,
    output:
        prcurve="{outdir}/results/model/train_test/{ref}_{db}/{model_id}/prcurve.pickle",
        plot="{outdir}/results/model/train_test/{ref}_{db}/{model_id}/prcurve.png",
    conda:
        "../envs/model.yml"
    log:
        "{outdir}/results/model/train_test/{ref}_{db}/{model_id}/prcurve.log",
    script:
        "../scripts/prcurve.py"
