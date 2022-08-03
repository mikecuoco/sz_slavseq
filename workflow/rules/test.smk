rule test_ref:
    input:
        expand("resources/{ref}{ext}", ext=extensions, ref=config["ref"])

rule test_cutadapt2:
    input:
        expand(["results/cutadapt2/{donor}/{type}/{sample}_R1.fastq.gz", "results/cutadapt2/{donor}/{type}/{sample}_R2.fastq.gz", "results/cutadapt2/{donor}/{type}/{sample}_qc.txt"],
                zip,
               sample=samples['sample_id'],
               donor=samples['donor_id'],
               type=samples['sample_type'])
               
rule test_bwa_mem:
    input:
        expand("results/bwa_mam/{donor}/{type}/{sample}.bam",
               zip,
               sample=samples['sample_id'],
               donor=samples['donor_id'],
               type=samples['sample_type'])

rule test_rmdup:
    input:
        expand("results/rmdup/{donor}/{type}/{sample}.bam",
               zip,
               sample=samples['sample_id'],
               donor=samples['donor_id'],
               type=samples['sample_type'])

rule test_tags:
    input:
        expand("results/tags/{donor}/{type}/{sample}.bam",
               zip,
               sample=samples['sample_id'],
               donor=samples['donor_id'],
               type=samples['sample_type'])

