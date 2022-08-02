rule test_ref:
    input:
        expand("resources/{ref}{ext}", ext=extensions, ref=config["ref"])

rule test_cutadapt2:
    input:
        expand(["results/cutadapt2/{sample}/{donor}_{type}_R1.fastq.gz", "results/cutadapt2/{sample}/{donor}_{type}_R2.fastq.gz", "results/cutadapt2/{sample}/{donor}_{type}_qc.txt"],
                zip,
               sample=samples['sample_id'],
               donor=samples['donor_id'],
               type=samples['sample_type'])
               
rule test_bwa_mem:
    input:
        expand("results/bwa_mam/{sample}/{donor}_{type}.bam",
               zip,
               sample=samples['sample_id'],
               donor=samples['donor_id'],
               type=samples['sample_type'])

rule test_rmdup:
    input:
        expand("results/rmdup/{sample}/{donor}_{type}.bam",
               zip,
               sample=samples['sample_id'],
               donor=samples['donor_id'],
               type=samples['sample_type'])

rule test_tags:
    input:
        expand("results/tags/{sample}/{donor}_{type}.bam",
               zip,
               sample=samples['sample_id'],
               donor=samples['donor_id'],
               type=samples['sample_type'])

