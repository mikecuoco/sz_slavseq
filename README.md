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
GENOME="hs37d5" # can be hs37d5, hs38dH or chm13v2
snakemake \
   all \
   --configfile .test/chr22/${GENOME}.yml \
   --cores 2 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --all-temp
```