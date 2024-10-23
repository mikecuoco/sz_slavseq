rule tune:
    input:
        data=rules.peaks_report.output,
    output:
        model=rules.peaks_report.output[0].replace(
            "data.bed.gz", "{features}/model.pkl"
        ),
        best_hp=rules.peaks_report.output[0].replace(
            "data.bed.gz", "{features}/best_hp.json"
        ),
        history=rules.peaks_report.output[0].replace(
            "data.bed.gz", "{features}/history.log"
        ),
        predictions=rules.peaks_report.output[0].replace(
            "data.bed.gz", "{features}/predictions.pqt"
        ),
        shap_values=rules.peaks_report.output[0].replace(
            "data.bed.gz", "{features}/shap_values.pqt"
        ),
    conda:
        "../envs/model_gpu.yml" if shutil.which("nvidia-smi") else "../envs/model.yml"
    params:
        max_iter=100,
        random_state=1,
        features=lambda wc: feature_sets[wc.features],
    threads: 1
    log:
        notebook=rules.peaks_report.log.notebook.replace(
            "peaks_report", "{features}/tune"
        ),
    notebook:
        "../scripts/tune.py.ipynb"


rule model_report:
    input:
        data=rules.regions_report.output,
        bulk=rules.germline_report.output,
        history=rules.tune.output.history,
        best_hp=rules.tune.output.best_hp,
        model=rules.tune.output.model,
    output:
        rules.tune.log.notebook.replace("tune", "model_report"),
    conda:
        "../envs/model.yml"
    params:
        random_state=1,
    log:
        notebook=rules.tune.log.notebook.replace("tune", "model_report"),
    notebook:
        "../scripts/model_report.py.ipynb"


rule calls_report:
    input:
        data=rules.tune.output.predictions,
    output:
        rules.model_report.log.notebook.replace("model_report", "calls_report"),
    conda:
        "../envs/model.yml"
    log:
        notebook=rules.model_report.log.notebook.replace("model_report", "calls_report"),
    notebook:
        "../scripts/calls_report.py.ipynb"


def get_donor_bams(wildcards):
    return expand(
        rules.sort.output,
        sample=samples.loc[wildcards.donor, "sample_id"],
        donor=wildcards.donor,
        allow_missing=True,
    )


rule make_igv_batch_script:
    input:
        regions=rules.tune.output.predictions.replace(
            "predictions.pqt", "{donor}.calls.bed"
        ),
        bams=get_donor_bams,
        rmsk=rules.run_rmsk.output.bed,
        megane_gaussian=expand(
            rules.vcf2bed.output, vcf="megane_gaussian", allow_missing=True
        ),
        primer_sites=rules.blast_primers.output.bed,
        fasta=config["genome"]["fasta"],
        fai=config["genome"]["fai"],
    output:
        rules.tune.output.predictions.replace(
            "predictions.pqt", "{donor}.igv_snapshots.bat"
        ),
    log:
        rules.tune.output.predictions.replace(
            "predictions.pqt", "{donor}.igv_snapshots.log"
        ),
    conda:
        "../envs/model.yml"
    params:
        maxPanelHeight=200,
        colorBy="READ_STRAND",
    script:
        "../scripts/igv_snapshots.py"


rule igv_download:
    output:
        "resources/IGV_Linux_2.17.4/igv_hidpi.sh",
    params:
        runtime="600",
        memory="1G",
    shell:
        """
        cd resources
        curl https://data.broadinstitute.org/igv/projects/downloads/2.17/IGV_Linux_2.17.4_WithJava.zip > IGV.zip 2> /dev/null
        unzip -o IGV.zip 2>&1 > /dev/null
        """


rule igv_snapshots:
    input:
        script=rules.make_igv_batch_script.output,
        igv=rules.igv_download.output,
    output:
        directory("{outdir}/results/{genome}/{reads}/peaks/{donor}/snapshots"),
    log:
        "{outdir}/results/{genome}/{reads}/peaks/{donor}/igv_snapshots.log",
    shell:
        """
        # requires xvfb-run to be installed by root for headless display
        xvfb-run -s "-screen 0, 7680x4320x24" --auto-servernum {input.igv} -b {input.script} &> {log}
        """
