#!/usr/bin/env python
# Created on: Jun 18, 2024 at 1:49:45 PM
__author__ = "Michael Cuoco"

from pathlib import Path
import pandas as pd

# example config
# config = {
# 	"genome": {
# 		"url": "https://example.com/genome.fa",
# 	},
# 	"tracks": [
# 		{
# 			"label": "track1",
# 			"url": "https://example.com/track1.bed",
# 		},
# 		{
# 			"label": "track2",
# 			"url": "https://example.com/track2.bed",
# 		},
# 		{
# 			"label": "track3",
# 			"url": "https://example.com/track3.bed",
# 		},
# 	]
# }

IGV = "/iblm/netapp/data3/mcuoco/sz_slavseq/resources/IGV_Linux_2.17.4/igv_hidpi.sh"
assert Path(IGV).exists(), f"{IGV} doesn't exists"

# TODO: make a class?
# TODO: make font size larger
def make_igv_batch_script(
    config: dict, regions: pd.DataFrame, outdir: str, flank: int = 5000
):
    """
    Generate a batch script for IGV to load a genome and tracks, and take snapshots of regions. Save images to outdir
    :param config: dictionary of URLs to track files
    :param regions: DataFrame with columns "Chromosome", "Start", "End"
    :param outdir: Directory to save snapshots
    :param flank_size: Number of bases to flank each region
    """

    # check inputs
    for c in ["Chromosome", "Start", "End"]:
        assert c in regions.columns, f"Column '{c}' not in regions DataFrame"
    assert "genome" in config.keys(), "Genome not found in config"
    assert "url" in config["genome"].keys(), "Genome URL not found in config"
    assert "tracks" in config.keys(), "Track URLs not found in config"
    for track in config["tracks"]:
        assert "label" in track.keys(), "Track label not found in config"
        assert "url" in track.keys(), "Track URL not found in config"

    genome = config["genome"]["url"]

    # make script
    script = outdir + "/igv.bat"
    with open(script, "w") as b:
        b.write(f"new\ngenome {genome}\nsnapshotDirectory {outdir}\n")  # type: ignore
        for track in config["tracks"]:
            b.write(f"load {track['url']} name={track['label']}\n")
            if track["url"].endswith(".bed"):
                b.write(f"expand {track['label']}\n")
                b.write(f"setTrackHeight {track['label']} 50\n")
            else:
                b.write(f"squish {track['label']}\n")
                b.write(f"setTrackHeight {track['label']} 200\n")

        # snapshot each row in regions
        for _, row in regions.iterrows():  # type: ignore
            start = row["Start"] - flank
            end = row["End"] + flank

            # write to batch script
            locus = f"{row['Chromosome']}:{start}-{end}"
            b.write(f"goto {locus}\n")
            b.write(f"snapshot {locus}.png\n")

        b.write(f"saveSession {outdir}/session.xml\n")
        b.write("exit\n")

    return script
