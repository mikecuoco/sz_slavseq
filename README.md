# SLAV-seq Snakemake Pipeline

[![Tests](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml/badge.svg)](https://github.com/mikecuoco/sz_slavseq/actions/workflows/main.yml)

Adapted from work by Apua Paquola and Ricardo Jacomini

## Development tips

If running on OSX, add `clang_osx-64` to `workflow/envs/align.yml` and ensure your `CONDA_BUILD_SYSROOT` environment variable points to a MACOSX sdk. See [conda-build compiler tools](https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html). This is necessary to install the [gapafim](https://github.com/apuapaquola/gapafim) dependency.

Below is the general development workflow

```bash
# format code
snakefmt .

# testing with specified reference genome and database
GENOME="hs37d5" # can be hs37d5, hs38dH or chm13v2
DB="eul1db" # can be eul1db or eul2db
snakemake \
   all \
   --configfile .test/chr21chr22/${GENOME}_${DB}.yml \
   --cores 2 \
   --use-conda \
   --show-failed-logs \
   --conda-cleanup-pkgs cache \
   --all-temp
```
