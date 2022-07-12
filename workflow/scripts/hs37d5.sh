#!/bin/bash

wget --read-timeout=300 -c --no-config -P "resources/" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz

samtools faidx hs37d5.fa

bwa index hs37d5.fa

cat resources/hs37d5.fa.fai | cut -f 1,2 > resources/hs37d5.genome