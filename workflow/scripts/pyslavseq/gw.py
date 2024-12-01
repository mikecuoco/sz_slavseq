#!/usr/bin/env python
# Created on: Nov 6, 2024 at 2:14:38 PM
__author__ = "Michael Cuoco"

# configure logging
import logging

logger = logging.getLogger(__name__)

from pathlib import Path
import subprocess
from tempfile import NamedTemporaryFile
import pandas as pd

RESULTS = Path(
    "/iblm/netapp/data3/mcuoco/sz_slavseq/results/chm13v2.0.XY/filtered/align"
)
RESOURCES = Path("/iblm/netapp/data3/mcuoco/sz_slavseq/resources/chm13v2.0.XY")
donors = (
    pd.read_csv("/iblm/netapp/data3/mcuoco/sz_slavseq/config/all_donors.tsv", sep="\t")[
        ["libd_id", "donor_id"]
    ]
    .set_index("libd_id")["donor_id"]
    .to_dict()
)

slavseq = {
    d: list((RESULTS / str(d)).rglob("*.tagged.sorted.bam")) for d in donors.values()
}
split_disc = {
    d: list((RESOURCES / f"wgs_calls/30x/{l}/").rglob("*split.disc.bw"))[0]
    for l, d in donors.items()
}
megane = {
    d: list((RESOURCES / f"wgs_calls/30x/{l}/").rglob("*MEI_final_gaussian.bed"))
    for l, d in donors.items()
}
l1hs = RESOURCES / "chm13v2.0.XY.fasta.rmsk.l1hs.bed"
ref = RESOURCES / "chm13v2.0.XY.fasta"


def run_gw(donor: str, regions: list, outdir):
    """
    Run GW on a list of BAM files
    :param bams: List of BAM files
    :param regions: List of regions to run GW on
    :param outdir: Output directory
    :param threads: Number of threads to use
    """

    # create tempfile
    with NamedTemporaryFile(suffix=".bed") as tmp:

        # make bedfile
        for r in regions:
            chrom, pos = r.split(":")
            start, end = pos.split("-")
            tmp.write(f"{chrom}\t{start}\t{end}\n".encode())

        # build command
        cmd = f"gw --no-show --theme igv"
        for b in slavseq[donor]:
            cmd += f" --bam {b}"
        cmd += f" --variants {tmp.name} --track {l1hs} --track {megane[donor]} --track {split_disc[donor]} --outdir {outdir} {ref}"

        # run
        logger.info(f"Running GW for {donor}")
        subprocess.run(cmd, shell=True, check=True)

    pass
