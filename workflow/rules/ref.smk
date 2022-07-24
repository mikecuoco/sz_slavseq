rule install_bwakit:
    output: directory("resources/bwa.kit")
    conda: "../envs/env.yml"
    shell:
        '''
        mkdir -p resources && cd resources
        wget -O- -q --no-config https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2 | tar xfj -
        '''

rule gen_ref:
    input: rules.install_bwakit.output
    output: "resources/{ref}/genome_og.fa"
    log: "resources/{ref}/gen_ref.log"
    conda: "../envs/env.yml"
    shell:
        '''
        touch {log} && exec 1>{log} 2>&1
        
        # run bwa.kit function
        {input}/run-gen-ref {wildcards.ref}
        mv {wildcards.ref}.fa {output}
        '''

rule fix_names_clean:
    input: 
        fa = rules.gen_ref.output,
        chrom_map = "workflow/scripts/hs37d5_to_hg19.tsv"
    output: 
        fa = "resources/{ref}/genome.fa",
        fai = "resources/{ref}/genome.fa.fai",
        chromsizes = "resources/{ref}/genome.genome"
    log: "resources/{ref}/fix_names.log"
    conda: "../envs/env.yml"
    shell:
        '''
        if [[ "{wildcards.ref}" =~ .*"37".*  ]]; then 
            # make hg19 names
            while read line; do
                hs37d5_chr=$(echo $line | cut -f1)
                hg19_chr=$(echo $line | cut -f2)
                sed -i "s/^>$hs37d5_chr/^>$hg19_chr/g" {input.fa} > {output.fa}
            done < {input.chrom_map}
        else
            mv {input.fa} {output.fa}
        fi

        samtools faidx {output.fa} 
        cut -f 1,2 {output.fai} > {output.chromsizes}
        '''

while read line; do
    hs37d5_chr=$(echo $line | cut -f1)
    hg19_chr=$(echo $line | cut -f2)
    sed -i "s/^>$hs37d5_chr/^>$hg19_chr/g" resources/hs37d5/genome_og.fa > resources/hs37d5/genome.fa
done < workflow/scripts/hs37d5_to_hg19.tsv

rule get_eul1db:
    input:
        chromsizes = expand("resources/{ref}/genome.genome",  ref=config["ref"])
    output:
        srip = "resources/eul1db/SRIP.txt",
        windows = "resources/eul1db/windows.csv"
    log: "resources/eul1db/eul1db.log"
    conda: "../envs/env.yml"
    shell:
        '''
        wget --no-config -q -O- http://eul1db.unice.fr/UserLists/DATA/downloads/SRIP.txt > {output.srip}
        workflow/scripts/eul1db_windows.py \
            --srip-file {output.srip} \
            --chromsizes {input.chromsizes} \
            --output {output.windows}
        '''

rule get_rmsk:
    output: "resources/{ref}/rmsk.txt.gz"
    log: "resources/{ref}/rmsk.log"
    conda: "../envs/env.yml"
    shell:
        '''
        touch {log} && exec 1>{log} 2>&1

        if [[ "{wildcards.ref}" =~ .*"38".*  ]]; then 
            url="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
        elif [[ "{wildcards.ref}" =~ .*"37".*  ]]; then 
            url="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz"
        fi

        wget --no-config -q -O- ${{url}} > {output}
        '''