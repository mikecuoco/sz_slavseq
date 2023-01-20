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
snakemake \
   all \
   --configfile .test/chr21chr22/${GENOME}.yml \
   --cores 2 \
   --use-conda \
   --show-failed-logs \
   --all-temp

# testing with SLURM on SSRDE
# NOTE: must increase --latency-wait when using a job scheduler

GENOME="hs37d5" # can be hs37d5, hs38dH or chm13v2
snakemake \
   all \
   --configfile .test/chr21chr22/${GENOME}.yml \
   --slurm \
   --jobs 8 \
   --default-resources slurm_partition=general runtime=60 \
   --resources cpus=50 mem_mb=10000 \
   --latency-wait 30 \
   --use-conda \
   --show-failed-logs
```

## Running >25 samples

NOTE: Ensure that storage device can handle highly parallel I/O. Set `--jobs` to be less than the number of cores on the storage device.

```bash
snakemake \
   all \
   --slurm \
   --jobs 300 \
   --default-resources slurm_partition=general \
   --resources cpus=1000 mem_mb=100000 \
   --latency-wait 30 \
   --use-conda \
   --show-failed-logs
```
