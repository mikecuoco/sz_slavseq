rule get_ref:
    output:
        multiext("resources/{ref}/genome", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".genome")
    log: "resources/{ref}/genome.log"
    conda: "../envs/env.yml"
    shell:
        '''
        touch {log} && exec 1>{log} 2>&1

        if [[ {wildcards.ref} == 'GRCh38' ]]; then 
            url="ftp://ftp.1000genomes.ebi.ac.uk:21/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla"
            extensions=(".fa" ".fa.amb" ".fa.ann" ".fa.bwt" ".fa.fai" ".fa.pac" ".fa.sa")

            # Download the reference genome, samtools index, and bwa index
            # If no data is received for more than 300 seconds, tell wget to resume the download
            for ext in ${{extensions[@]}}; do
                while true; do
                    wget -q --read-timeout=300 -c --no-config -P "resources/{wildcards.ref}" -O "resources/{wildcards.ref}/genome${{ext}}" "${{url}}${{ext}}" && break
                done
            done

        elif [[ {wildcards.ref} == 'hs37d5' ]]; then
            url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"

            # Download the reference genome only, then generate indices ourselves
            # If no data is received for more than 300 seconds during download, tell wget to resume the download
            wget -q --read-timeout=300 -c --no-config -P "resources/{wildcards.ref}" -O "resources/{wildcards.ref}/genome${{ext}}" {{url}}
            gunzip resources/{wildcards.ref}/genome.fa.gz
            
            samtools faidx resources/{wildcards.ref}/genome.fa
            bwa index resouces/{wildcards.ref}/genome.fa
        fi

        cat resources/{wildcards.ref}/genome.fa.fai | cut -f 1,2 > resources/{wildcards.ref}/genome.genome
        '''

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
        wget --no-config -q -P resources/eul1db/ http://eul1db.unice.fr/UserLists/DATA/downloads/SRIP.txt
        workflow/scripts/eul1db_windows.py \
            --srip-file {output.srip} \
            --chromsizes {input.chromsizes} \
            --output {output.windows}
        '''