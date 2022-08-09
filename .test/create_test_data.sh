#!/bin/bash
#
# Author: Rohini Gadde
# 
# Use this script to create new test data.
# 
# Usage: ./create_test_data.sh <bam> <fastq> <out_prefix> <region>
# bam = input bam file
# fastq = path to corresponding fastq file (either mate in pair)
# out_prefix = prefix for the output files (including outdir)
# region = region of the genome in format (chrN:start-end)

if [[ $# -lt 3 ]] || [[ $# -gt 4 ]]; then 
    echo "Usage: ./create_test_data.sh <bam> <fastq> <out_prefix> <region>." \
        "First three arguments necessary!" \
        "FASTQ argument is a path to a single file (can be either read in pair)"
    exit 1
fi

export LC_ALL=C

outdir=$(git rev-parse --show-toplevel)
outdir+="/.test/test_fastqs"

if [[ "$2" == *R1* ]]; then
    mate=$(echo "$2" | sed 's/R1/R2/')
elif [[ "$2" == *R2* ]]; then
    mate=$(echo "$2" | sed 's/R2/R1/')
else
    echo "FASTQ incorrectly named"
    exit 1
fi

echo "$2"
echo "$mate"

if [[ $# -eq 3 ]]; then
    samtools view "$1" | awk '{print $1}' > tmp_ids.txt
else
    samtools view "$1" "$4" | awk '{print $1}' > tmp_ids.txt
fi

zcat "$2" | fgrep -f tmp_ids.txt -A 3 --no-group-separator - > "$3_R1.fastq"
zcat "${mate}" | fgrep -f tmp_ids.txt -A 3 --no-group-separator - > "$3_R2.fastq"

rm tmp_ids.txt
