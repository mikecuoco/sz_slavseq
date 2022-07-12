#!/bin/bash

url="ftp://ftp.1000genomes.ebi.ac.uk:21/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla"
extensions=(".fa" ".fa.amb" ".fa.ann" ".fa.bwt" ".fa.fai" ".fa.pac" ".fa.sa")

# Download the reference genome, samtools index, and bwa index
# If no data is received for more than 300 seconds during download, tell wget to resume the download
for ext in ${extensions[@]}; do
    while true; do
        wget --read-timeout=300 -c --no-config -P resources/ "${url}${ext}" && break
    done
done

cat resources/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai | cut -f 1,2 > resources/hg38.genome
