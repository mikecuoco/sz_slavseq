#!/usr/bin/bash
#
# Author: Rohini Gadde
# Usage: ./liftover.sh
# 
# Liftover L1 (single-base) insertion positions from hg19 to hg38

if [[ "$6" == "hg19ToHg38.over.chain.gz" ]]; then
    echo "True"
    wget -P resources/ -q --no-config \
        https://raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh37.novel_CUPs.bed

    # bed -> vcf
    Rscript workflow/scripts/bed_to_vcf.R "$1" \
        "resources/$4/insertions.vcf" \
        "resources/$5/genome.fa"

    # Remove unstable positions
    vcftools --vcf insertions.vcf --exclude-bed resources/FASTA_BED.ALL_GRCh37.novel_CUPs.bed \
        --recode --recode-INFO-all --out insertions_stable.vcf
    
    # picard LiftoverVcf
    picard CreateSequenceDictionary -R "$2" -O "resources/$5/genome.dict"
    picard LiftoverVcf -I insertions_stable.vcf -O insertions_lifted.vcf -C "resources/$6" \
        --REJECT rejected_insertions.vcf -R "$2"

    # vcf -> bed
    tail -n +6 "resources/$4/insertions_lifted.vcf" | \
        awk -v OFS="\t" '{print $1, $2-1, $2}' \
        > liftover.tmp
    sed -i '1s/^/chr\tstart\tend\n/' liftover.tmp > "$3"
    rm liftover.tmp
else
    # no liftover necessary
    echo "False"
    mv "$1" "$3"
fi