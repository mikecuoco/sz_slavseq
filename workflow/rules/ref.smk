rule get_ref:
    output:
        multiext("resources/{ref}", ".fa", ".fa.amb", ".fa.ann", ".fa.bwt", ".fa.fai", ".fa.pac", ".fa.sa", ".genome")
    log:
        "resources/{ref}.log"
    conda:
        "../envs/env.yml"
    shell:
        '''
        touch {log} && exec 1>{log} 2>&1

        if [[ {wildcards.ref} == 'GRCh38' ]]; then 
            url="ftp://ftp.1000genomes.ebi.ac.uk:21/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla"
            extensions=(".fa" ".fa.amb" ".fa.ann" ".fa.bwt" ".fa.fai" ".fa.pac" ".fa.sa")

            # Download the reference genome, samtools index, and bwa index
            # If no data is received for more than 300 seconds during download, tell wget to resume the download
            for ext in ${{extensions[@]}}; do
                while true; do
                    wget -q --read-timeout=300 -c --no-config -P "resources/" -O "GRCh38${{ext}}" "${{url}}${{ext}}" && break
                done
            done

            cat resources/GRCh38.fa.fai | cut -f 1,2 > resources/GRCh38.genome
        elif [[ {wildcards.ref} == 'hs37d5' ]]; then
            # If no data is received for more than 300 seconds during download, tell wget to resume the download
            wget -q --read-timeout=300 -c --no-config -P "resources/" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
            gunzip resources/hs37d5.fa.gz

            samtools faidx resources/hs37d5.fa

            bwa index resources/hs37d5.fa

            cat resources/hs37d5.fa.fai | cut -f 1,2 > resources/hs37d5.genome
        fi
        '''