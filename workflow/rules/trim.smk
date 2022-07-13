rule cutadapt:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample, "R1"],
        r2 = lambda wildcards: samples.loc[wildcards.sample, "R2"]
    output:
        r1 = 'results/cutadapt/{sample}/{donor}_{type}_R1.fastq.gz',
        r2 = 'results/cutadapt/{sample}/{donor}_{type}_R2.fastq.gz'
    log:
        log1 = "results/cutadapt/{sample}/{donor}_{type}.log",
        log2 = "results/cutadapt/{sample}/{donor}_{type}.log"
    threads: 4
    conda: "../envs/cutadapt.yml"
    shell:
        ''' 
        function rc () {{
            echo $1 | tr ACGTacgt TGCAtgca | rev
        }}
        export -f rc

        R1_ADAPTER='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
        R2_ADAPTER='CAAGCAGAAGACGGCATACGAGANNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'

        NESTED_PRIMER='TAACTAACCTGCACAATGTGCAC'

        R1_FRONT=${{R1_ADAPTER}}
        R2_FRONT=${{R2_ADAPTER}}${{NESTED_PRIMER}}
        R1_END=`rc ${{R2_FRONT}}`
        R2_END=`rc ${{R1_FRONT}}`

        QUALITY_BASE=33
        QUALITY_CUTOFF=28
        MINIMUM_LENGTH=36
        ADAPTOR_OVERLAP=5
        ADAPTOR_TIMES=4
        
        cutadapt -j {threads} \
            --quality-base=${{QUALITY_BASE}} \
            --quality-cutoff=${{QUALITY_CUTOFF}} \
            --minimum-length=${{MINIMUM_LENGTH}} \
            --overlap=${{ADAPTOR_OVERLAP}} \
            --times=${{ADAPTOR_TIMES}} \
            --front=${{R1_FRONT}} \
            --adapter=${{R1_END}} \
            --paired-output tmp.2.{wildcards.sample}.fastq -o tmp.1.{wildcards.sample}.fastq {input.r1} {input.r2} \
            > {log.log1} 2>&1

        cutadapt -j {threads} \
            --quality-base=${{QUALITY_BASE}} \
            --quality-cutoff=${{QUALITY_CUTOFF}} \
            --minimum-length=${{MINIMUM_LENGTH}} \
            --overlap=${{ADAPTOR_OVERLAP}} \
            --times=${{ADAPTOR_TIMES}} \
            --front=${{R2_FRONT}} \
            --adapter=${{R2_END}} \
            --paired-output {output.r1} -o {output.r2} tmp.2.{wildcards.sample}.fastq tmp.1.{wildcards.sample}.fastq \
            > {log.log2} 2>&1

        rm -f tmp.2.{wildcards.sample}.fastq tmp.1.{wildcards.sample}.fastq
        '''    