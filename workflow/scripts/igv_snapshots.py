#!/usr/bin/env python
# Created on: Apr 15, 2024 at 8:08:03 PM
__author__ = "Michael Cuoco"

import logging

logging.basicConfig(
    filename=snakemake.log[0],  # type: ignore
    filemode="w",
    level=logging.INFO,
)

from pathlib import Path
import pyranges as pr
from snakemake.shell import shell

batch_script_fn = Path(snakemake.output[0]) / "igv_snapshots.bat"  # type: ignore
bamlist_fn = Path(snakemake.output[0]) / "bamlist.txt"  # type: ignore
batch_script_fn.parent.mkdir(parents=True, exist_ok=True)
with open(batch_script_fn, "w") as b, open(bamlist_fn, "w") as l:
    b.write(f"new\ngenome {snakemake.input.fasta}\nsnapshotDirectory {snakemake.output[0]}\n")  # type: ignore
    # b.write(f"new\ngenome hs1\nsnapshotDirectory {snakemake.output[0]}\n") # type: ignore

    b.write(f"load {snakemake.input.megane_gaussian}\n")  # type: ignore
    b.write(f"load {snakemake.input.rmsk}\n")  # type: ignore
    b.write(f"load {snakemake.input.primer_sites}\n")  # type: ignore

    for bam in snakemake.input.bams:  # type: ignore
        l.write(bam + "\n")

    b.write(f"load {str(bamlist_fn)}\n")
    b.write(f"maxPanelHeight {snakemake.params.maxPanelHeight}\n")  # type: ignore

    for _, row in pr.read_bed(snakemake.input.regions).df.iterrows():  # type: ignore
        locus = f"{row['Chromosome']}:{row['Start']-1000}-{row['End']+1000}"
        b.write(f"goto {locus}\n")
        b.write(f"colorBy {snakemake.params.colorBy}\n")  # type: ignore
        b.write(f"snapshot {locus}.png\n")

    b.write("exit\n")

shell(
    'xvfb-run --auto-servernum --server-num=1 java -Xmx4000m --module-path="${{CONDA_PREFIX}}/bin/igv" -b {batch_script_fn} &> {snakemake.log[0]}'
)
