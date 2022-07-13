rule test_ref:
    input:
        expand("resources/{ref}{ext}", ext=extensions, ref=config["ref"])
        
rule test_cutadapt:
    input:
        expand(["results/cutadapt/{sample}/{donor}_{type}_R1.fastq.gz","results/cutadapt/{sample}/{donor}_{type}_R2.fastq.gz"],
               zip,
               sample=samples['sample_id'],
               donor=samples['donor_id'],
               type=samples['sample_type'])