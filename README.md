# SLAV-seq Snakemake Pipeline

[![Tests](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml/badge.svg)](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml)

Adapted from work by Apua Paquola and Ricardo Jacomini

## Development tips

```bash
# format code
snakefmt .

# run lint checks
snakemake --lint

# Testing with hs37d5 reference (can also use configs for hs38dH and chm13v2)
# NOTE: Tests do not work from within tmux environment
# Set cachedir
export SNAKEMAKE_OUTPUT_CACHE=$(pwd)/.snakemake-cache
mkdir -p $SNAKEMAKE_OUTPUT_CACHE

# Quick test on steps through feature extraction
GENOME="hs37d5" # can be hs37d5, hs38dH or chm13v2
TARGET="results/folds/${GENOME}/test/mda/Y_test.pickle.gz"
snakemake \
   $TARGET \
   --configfile .test/chr22/${GENOME}.yml \
   --cores 2 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --all-temp

# Quick test on steps following making folds
mkdir -p results/folds/ && cp -r .test/data/${GENOME} results/folds/
TARGET="all"
snakemake \
   $TARGET \
   --configfile .test/allchrs/${GENOME}.yml \
   --cores 2 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --rerun-triggers mtime \
   --all-temp 
```