# SLAV-seq Snakemake Pipeline

![Tests](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml/badge.svg)

Adapted from work by Apua Paquola and Ricardo Jacomini

## Development tips

```bash
# format code
snakefmt .

# run lint checks
snakemake --lint

# Testing
# NOTE: Tests do not work from within tmux environment
# Set cachedir
export SNAKEMAKE_OUTPUT_CACHE=$(pwd)/.snakemake-cache
mkdir -p $SNAKEMAKE_OUTPUT_CACHE

# Quick test on steps through feature extraction
snakemake \
   --configfile .test/config_chr22/config.yml \
   --until results/flank_features/test/mda/chr22.pickle.gz \
   --cores 4 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --all-temp \
   --cache

# Quick test on steps following feature extraction
snakemake \
   --configfile .test/config_model/config.yml \
   all \
   --cores 4 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --all-temp \
   --cache

# Test the pipeline all the way through
snakemake \
   --configfile .test/config_sample/config.yml \
    all \
   --cores 4 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --all-temp \
   --cache
```