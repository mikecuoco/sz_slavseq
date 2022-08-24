# SLAV-seq Snakemake Pipeline

![Tests](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml/badge.svg)

Adapted from work by Apua Paquola and Ricardo Jacomini

## Development tips

```bash
# format code
snakefmt .

# run lint checks
snakemake --lint

# test a single sample
snakemake --cores 1 all --configfile .test/config_sample/config.yml --rerun-incomplete --show-failed-logs --use-conda --debug

# test multiple samples from a single individual
snakemake --cores 8 all --configfile .test/config_indv/config.yml --rerun-incomplete --show-failed-logs --use-conda 
```