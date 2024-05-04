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

batch_script_fn = snakemake.output[0]  # type: ignore
outdir = Path(batch_script_fn).parent
Path(outdir).mkdir(parents=True, exist_ok=True)
regions = pr.read_bed(snakemake.input.regions).merge().df  # type: ignore

with open(batch_script_fn, "w") as b:
    b.write(f"new\ngenome {snakemake.input.fasta}\nsnapshotDirectory {outdir}/snapshots\n")  # type: ignore

    hip, dlpfc = [], []
    for bam in snakemake.input.bams:  # type: ignore
        if "gDNA" in Path(bam).name:
            bulk_bam = bam
        elif "ush" in Path(bam).name.lower():
            hip.append(bam)
        elif "usd" in Path(bam).name.lower():
            dlpfc.append(bam)

    for bam in hip:
        name = Path(bam).name.split(".")[0]
        b.write(f"load {bam} name={name}\n")
    for bam in dlpfc:
        name = Path(bam).name.split(".")[0]
        b.write(f"load {bam} name={name}\n")

    b.write(f"maxPanelHeight {snakemake.params.maxPanelHeight}\n")  # type: ignore
    b.write(f"load {snakemake.input.megane_gaussian} name=wgs_calls\n")  # type: ignore
    b.write(f"load {snakemake.input.rmsk} name=RepeatMasker_LINE-1\n")  # type: ignore
    b.write(f"load {snakemake.input.primer_sites} name=predicted_SLAVseq_primer_sites\n")  # type: ignore
    b.write(f"load {snakemake.input.regions} name=SLAVseq_calls\n")  # type: ignore
    name = Path(bulk_bam).name.split(".")[0]
    b.write(f"load {bulk_bam} name={name}\n")  # type: ignore
    b.write("squish\n")

    for _, row in regions.iterrows():  # type: ignore
        locus = f"{row['Chromosome']}:{row['Start']-1000}-{row['End']+1000}"
        b.write(f"goto {locus}\n")
        b.write(f"colorBy {snakemake.params.colorBy}\n")  # type: ignore
        b.write(f"snapshot {locus}.png\n")

    b.write("exit\n")
