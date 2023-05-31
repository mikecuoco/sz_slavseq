# SLAV-seq Snakemake Pipeline

[![Tests](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml/badge.svg)](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml)

Adapted from work by Apua Paquola and Ricardo Jacomini

## Development tips

Below is the general development workflow

```bash
# install precommit hooks
pre-commit install

# test
snakemake \
   all \
   --configfile .test/chr21chr22/config.yml \
   --cores 2 \
   --use-conda \
   --show-failed-logs \
   --all-temp

# run the workflow
snakemake \
   all \
   --cores 8 \
   --use-conda \
   --show-failed-logs \
   --until run_rmsk
```

## TODO:

05/17/2023

- [ ] Get left and right bg windows in feature table in `./workflow/scripts/get_features.py`
- [ ] Keep RepeatMasker insertions only if position in end > 860
- [ ] Only consider read pairs with at least one mate having MAPQ = 60
- [ ] Call germline L1 from bulk SLAV-seq data, validate with WGS
- [ ] Remove mean template length, instead count # of discordant and concordant reads
- [ ] switch back to T2T-CHM13
- [ ] Try peak calling again
- [ ] organize notebooks into subfolders with python scripts for utility functions
