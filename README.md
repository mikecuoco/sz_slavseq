# SLAV-seq Snakemake Pipeline

[![Tests](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml/badge.svg)](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml)

Adapted from work by Apua Paquola and Ricardo Jacomini

## Development tips

Below is the general development workflow

```bash
# install precommit hooks
pre-commit install

# test, requires proper setup of .test/resources
snakemake \
   all \
   --directory .test \
   --cores 2 \
   --use-conda \
   --show-failed-logs \
   --all-temp

# run the workflow
snakemake \
   all \
   --cores 8 \
   --use-conda \
   --show-failed-logs
```
