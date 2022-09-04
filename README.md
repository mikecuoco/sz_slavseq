# SLAV-seq Snakemake Pipeline

![Tests](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml/badge.svg)

Adapted from work by Apua Paquola and Ricardo Jacomini

## Development tips

```bash
# format code
snakefmt .

# run lint checks
snakemake --lint

# Test
# NOTE: Tests do not work from within tmux environment
# Set cachedir
export SNAKEMAKE_OUTPUT_CACHE=$(pwd)/.snakemake-cache
mkdir -p $SNAKEMAKE_OUTPUT_CACHE
# Run test
snakemake --cores 1 all --configfile .test/config_sample/config.yml --cache --rerun-incomplete --show-failed-logs --use-conda --debug
```