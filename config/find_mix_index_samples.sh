#!/usr/bin/env bash
# Author: Mike Cuoco
# Created on: 9/6/24, 1:43 PM


# exit if any non-zero, exit if undefined var
set -euo pipefail

paths=(
	"/iblm/logglun01/fqiu/AWS/for-ndar/SLAV-Seq/"
	"/iblm/logglun02/BSMN/SLAV-Seq/"
	"/iblm/netapp/data3/mcuoco/sz_slavseq/data/"
)

if [ -f "indices.txt" ]; then
	rm -f indices.txt
fi

# extract flowcell and index from read
function uniq_indices {
	echo "Finding unique indices for $1"
	filename=$(basename "$1")
	inds=$(zcat "$1" | awk "NR%4==1" | awk '{split($1, a, ":"); print a[1],a[3],a[4],$2}' | sort -u | tr '\n' ',')
	printf "$1\t$inds\n" >> indices.txt
}
export -f uniq_indices

# get all the files
files=()
for p in "${paths[@]}"; do
	files+=($(find "$p" -name "*.fastq.gz"))
done

echo "Found ${#files[@]} fastq files"

parallel --progress -j 8 'uniq_indices {}' ::: "${files[@]}"
