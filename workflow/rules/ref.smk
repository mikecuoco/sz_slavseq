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

# TODO: edit fix_names.py to also change the names for hg38
rule fix_names_clean:
    input: 
        fa = rules.gen_ref.output
    output: 
        fa = "resources/{ref}/genome.fa",
        fai = "resources/{ref}/genome.fa.fai",
        chromsizes = "resources/{ref}/genome.genome"
    log: "resources/{ref}/fix_names.log"
    conda: "../envs/env.yml"
    script: "../scripts/fix_names.py"
        
rule get_eul1db:
    input: expand("resources/{ref}/genome.genome",  ref=config["ref"])
    output: "resources/eul1db/windows.csv"
    conda: "../envs/env.yml"
    script: "../scripts/get_eul1db.py"