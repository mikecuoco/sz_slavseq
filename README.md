# SLAV-seq Snakemake Pipeline

![Tests](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml/badge.svg)

Adapted from work by Apua Paquola and Ricardo Jacomini

## Development tips

Before pushing on the main branch run the following from the project root:

```bash
snakemake --lint
snakemake --cores 1 all --configfile .test/config/config.yml --rerun-incomplete --show-failed-logs --use-conda --debug
```
