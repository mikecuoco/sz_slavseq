#!/bin/bash
#
# Author: Rohini Gadde
# 
# Use this script to create new test data.
# 
# Usage: ./create_test_data.sh <bam> <fastq_R1> <fastq_R2> <out_prefix> <region>
# bam = input bam file
# fastq_R1 = path to corresponding R1 fastq file
# fastq_R2 = path to corresponding R2 fastq file
# out_prefix = prefix for the output files (including outdir)
# region = region of the genome in format (chrN:start-end)

export LC_ALL=C

if [[ $# -lt 4 ]]; then 
    echo "Usage: ./create_test_data.sh <bam> <fastq_R1> <fastq_R2> <out_prefix> <region1> [<region2> ...]." \
        "First four arguments necessary!"
    exit 1
fi

# put regions into an array
if [[ $# -gt 4 ]]; then
    regions=()
    for n in $(seq 5 $#); do
        regions+=("${!n}")
    done
fi

if [[ "$2" != *R1*.fastq* ]] || [[ "$3" != *R2*.fastq* ]]; then
    echo "FASTQs incorrectly named"
    exit 1
fi

if [[ $# -eq 4 ]]; then
    samtools view "$1" | awk '{print $1}' > tmp_ids.txt
else
    samtools view "$1" "${regions[@]}" | awk '{print $1}' > tmp_ids.txt
fi

zcat "$2" | fgrep -f tmp_ids.txt -A 3 --no-group-separator - > "$4_R1.fastq"
zcat "$3" | fgrep -f tmp_ids.txt -A 3 --no-group-separator - > "$4_R2.fastq"
gzip $4_R*.fastq

rm tmp_ids.txt
