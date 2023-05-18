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
- [ ] Write read2 in `./workflow/scripts/tag_reads.py` if it uniquely maps or is in proper pair
- [ ] Only consider reads with at least one mate having MAPQ = 60
- [ ] Call germline L1 from bulk SLAV-seq data, validate with WGS
- [ ] Remove mean template length, instead count # of discordant and concordant reads
