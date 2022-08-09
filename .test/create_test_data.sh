#!/bin/bash
#
# Author: Rohini Gadde
# 
# Use this script to create new test data.
# 
# Usage: ./create_test_data.sh <bam> <fastq> <out_prefix> <region>
# bam = input bam file
# fastq_R1 = path to corresponding R1 fastq file
# fastq_R2 = path to corresponding R2 fastq file
# out_prefix = prefix for the output files (including outdir)
# region = region of the genome in format (chrN:start-end)

if [[ $# -lt 4 ]] || [[ $# -gt 5 ]]; then 
    echo "Usage: ./create_test_data.sh <bam> <fastq_R1> <fastq_R2> <out_prefix> <region>." \
        "First four arguments necessary!"
    exit 1
fi

export LC_ALL=C

if [[ "$2" != *R1.fastq* ]] || [[ "$3" != *R2.fastq* ]]; then
    echo "FASTQs incorrectly named"
    exit 1
fi

if [[ $# -eq 4 ]]; then
    samtools view "$1" | awk '{print $1}' > tmp_ids.txt
else
    samtools view "$1" "$5" | awk '{print $1}' > tmp_ids.txt
fi

zcat "$2" | fgrep -f tmp_ids.txt -A 3 --no-group-separator - > "$4_R1.fastq"
zcat "$3" | fgrep -f tmp_ids.txt -A 3 --no-group-separator - > "$4_R2.fastq"

rm tmp_ids.txt
