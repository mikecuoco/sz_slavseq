rule get_ref:
    output: "resources/{ref}/genome.fa"
    log: "resources/{ref}/get_ref.log"
    conda: "../envs/env.yml"
    shell:
        '''
        touch {log} && exec 1>{log} 2>&1

        if [[ {wildcards.ref} == 'GRCh38' ]]; then 
            url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
        elif [[ {wildcards.ref} == 'hs37d5' ]]; then
            url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
        fi

        # If no data is received for more than 300 seconds during download, tell wget to resume the download
        wget -q --read-timeout=300 -c --no-config -P "resources/{wildcards.ref}" -O "{output}.gz" ${{url}}
        gunzip {output}.gz
        '''

# TODO: perform liftover depending on ref
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

rule get_rmsk:
    output: "resources/{ref}/rmsk.txt.gz"
    log: "resources/{ref}/rmsk.log"
    conda: "../envs/env.yml"
    shell:
        '''
        touch {log} && exec 1>{log} 2>&1

        if [[ {wildcards.ref} == 'GRCh38' ]]; then 
            url="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
            wget --no-config -q -P resources/{wildcards.ref} -O resources/{wildcards.ref}/rmsk.txt.gz ${{url}}
        elif [[ {wildcards.ref} == 'hs37d5' ]]; then
            url="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz"
            wget --no-config -q -P resources/{wildcards.ref} -O resources/{wildcards.ref}/rmsk.txt.gz ${{url}}
        fi
        '''