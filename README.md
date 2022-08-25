# SLAV-seq Snakemake Pipeline

![Tests](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml/badge.svg)

Adapted from work by Apua Paquola and Ricardo Jacomini

## Development tips

```bash
# format code
snakefmt .

# run lint checks
snakemake --lint

# NOTE: sra download for tests fails within tmux session
# test a single sample
export SNAKEMAKE_OUTPUT_CACHE=$(pwd)/.snakemake-cache
snakemake --cores 1 all --configfile .test/config_sample/config.yml --rerun-incomplete --show-failed-logs --use-conda --debug

# test multiple samples from a single individual
snakemake --cores 8 all --configfile .test/config_indv/config.yml --rerun-incomplete --show-failed-logs --use-conda 
```