#!/usr/bin/env bash
# Author: Mike Cuoco
# Liftover input bed file from one genome build to another using CrossMap.py in snakemake

# exit if any non-zero, exit if undefined var
set -euo pipefail

# set input vars and logging
touch ${snakemake_log} && exec 2>${snakemake_log} 
source=${snakemake_params[source]}
target=${snakemake_params[target]}
input=${snakemake_input[0]}
output=${snakemake_output[0]}

# select chain file based off input source and target
# TODO: look for CUP Files for other liftovers
# TODO: can we use source CUP files for liftovers to other targets?
if [[ "$source" == "hg19" ]]; then

	if [[ "$target" == *"38"* ]]; then
		CUP_URL="https://raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh37.novel_CUPs.bed"
		URL="https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz"
	elif [[ "$target" == "chm13v2" ]]; then
		URL="https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg19-chm13v2.over.chain.gz"
	fi

elif [[ "$source" == *"38"* ]]; then

	if [[ "$target" == "hg19" ]]; then
		CUP_URL="https://raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh38.novel_CUPs.bed"
		URL="https://hgdownload.soe.ucsc.edu/gbdb/hg38/liftOver/hg38ToHg19.over.chain.gz"
	elif [[ "$target" == "chm13v2" ]]; then
		CUP_URL="https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38.liftover-mask.bed"
		URL="https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/hg38-chm13v2.over.chain.gz"
	fi
	
elif [[ "$source" == "chm13v2" ]]; then

	if [[ "$target" == "hg19" ]]; then
		URL="https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/chm13v2-hg19.over.chain.gz"
	elif [[ "$target" == *"38"* ]]; then
		URL="https://hgdownload.gi.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/liftOver/chm13v2-hg38.over.chain.gz"
	fi

else

	echo "unable to liftover, chain file unavailable" && exit 1

fi

# download Conversion Unstable Positions (CUPs) to exclude from liftover input
final_input=$input
# if [[ -v CUP_URL ]]; then
# 	wget -O resources/${source}To${target}_mask.bed -q --no-config $CUP_URL

# 	final_input=$(mktemp)
# 	grep -v -f resources/${source}To${target}_mask.bed $input > $final_input
# fi

# download chain file
CHAINFILE="resources/$(basename $URL)"
mkdir -p $(dirname $CHAINFILE) # make dir if doesn't exist
wget -O $CHAINFILE -q --no-config $URL

# run liftOver using CrossMap wrapper
CrossMap.py bed $CHAINFILE $final_input $output

exit 0
