#!/usr/bin/env bash
# Download the reference genome and decompress 
# Usage: bash run-gen-ref.sh <chm13v2|hs38|hs38a|hs38DH|hs37|hs37d5>
# TODO: find bwa indices and download them

# store urls for downloads, 13v2 and 38 are from NCBI GenBank
url13v2="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz"
# url38="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"
url38="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz"
url37d5="ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
urlbwakit="https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2"

if [ $# -eq 0 ]; then
	echo "Usage: $0 <chm13v2|hs38|hs38a|hs38DH|hs37|hs37d5>"
	echo "Analysis sets:"
	echo "  chm13v2  complete genome from CHM13htert cell line"
	echo "  hs38     primary assembly of GRCh38 (incl. chromosomes, unplaced and unlocalized contigs) and EBV"
	echo "  hs38a    hs38 plus ALT contigs"
	echo "  hs38DH   hs38a plus decoy contigs and HLA genes (recommended for GRCh38 mapping)"
	echo "  hs37     primary assembly of GRCh37 (used by 1000g phase 1) plus the EBV genome"
	echo "  hs37d5   hs37 plus decoy contigs (used by 1000g phase 3)"
	echo ""
	exit 1;
fi

# if download hs38, make temp dir and download bwa.kit
if [[ "$1" == *"hs38"* ]]; then
	tmp=$(mktemp -d)
	cd $tmp
	trap 'rm -rf -- "$tmp"' EXIT

	wget -O- $urlbwakit | tar xfj -
	cd -
fi

if [ $1 == "chm13v2" ]; then
	wget -O- $url13v2 | gzip -dc > $1.fa 2> /dev/null
elif [ $1 == "hs38DH" ]; then
	(wget -O- $url38 | gzip -dc; cat $tmp/bwa.kit/resource-GRCh38/hs38DH-extra.fa) > $1.fa 
	# [ ! -f $1.fa.alt ] && cp $tmp/bwa.kit/resource-GRCh38/hs38DH.fa.alt $1.fa.alt
elif [ $1 == "hs38a" ]; then
	wget -O- $url38 | gzip -dc > $1.fa
	# [ ! -f $1.fa.alt ] && grep _alt $tmp/bwa.kit/resource-GRCh38/hs38DH.fa.alt > $1.fa.alt
elif [ $1 == "hs38" ]; then
	wget -O- $url38 | gzip -dc | awk '/^>/{f=/_alt/?0:1}f' > $1.fa
elif [ $1 == "hs37d5" ]; then
	wget -O- $url37d5 | gzip -dc > $1.fa 2>/dev/null
elif [ $1 == "hs37" ]; then
	wget -O- $url37d5 | gzip -dc 2>/dev/null | awk '/^>/{f=/>hs37d5/?0:1}f' > $1.fa
else
	echo "ERROR: unknown genome build"
fi

