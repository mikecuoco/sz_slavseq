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
   --configfile .test/config/config.yml \
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

- [ ] switch back to T2T-CHM13
- [ ] switch to parquet files for features and labels outputs
- [ ] Correct for GC bias?

Filter out regions (see RetroSom paper)

- known insertions +/- 1kb
- low complexity regions
- GRC blacklist regions
- 10x blacklist regions
- decoy regions

Features to add:

- separate classes for transductions and regular insertions? (transductions have discordant read2)
- start/end sites of clipped reads
- alignment score quanties as fraction of read length
- number of supplementary alignments / read
- fraction of proper pairs
- L1 consensus start site
- L1 consensus AS
- ORF2 motif
