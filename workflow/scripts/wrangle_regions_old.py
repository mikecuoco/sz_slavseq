#!/usr/bin/env python
# Created on: Oct 4, 2024 at 9:58:54 AM
__author__ = "Michael Cuoco"

import logging

logger = logging.getLogger(__name__)
from time import time
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
import pickle as pkl

import pandas as pd
import numpy as np
import pyranges as pr
from scipy.sparse import csr_matrix
import anndata as ad

from Bio import motifs, SeqIO
from Bio.Seq import Seq
from pyfaidx import Fasta
from pyslavseq.preprocessing import collate_labels


def compute_distance(df: pd.DataFrame, reference_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute distance of each region to the nearest germline region
    df: pd.DataFrame, intervals
    germline_df: pd.DataFrame, germline regions
    """
    df = pr.PyRanges(df).sort().df
    gdf = pr.PyRanges(germline_df)
    donor_df = (
        pr.PyRanges(df)
        .nearest(gdf, overlap=False)
        .df.rename(columns={"Distance": "germline_distance"})
    )
    donor_df["germline_distance"] = donor_df["germline_distance"].abs()
    donor_df.drop(columns=["Start_b", "End_b"], inplace=True)
    return donor_df


def merge_clusters(cluster: int, df: pd.DataFrame, labels: list) -> pd.Series:
    """
    Merge clusters and compute motif scores, groupby apply function
    :param cluster: int, cluster number
    :param df: pd.DataFrame, cluster regions
    :param labels: list, labels
    """

    # get coordinates and sequence of this cluster
    chrom = df["Chromosome"].iloc[0]
    start = df["Start"].min()
    end = df["End"].max()
    return pd.Series(
        {
            "Cluster": cluster,
            "locus": f"{chrom}:{start}-{end}",
            "Chromosome": chrom,
            "Start": start,
            "End": end,
            "width": end - start,
            "n_ref_reads": df["n_ref_reads"].max(),
            "max_mapq": df["mapq_q1"].max(),
            "n_reads": df["n_reads"].sample(random_state=42).iloc[0],
            "three_end_clippedA": df["three_end_clippedA_q1"]
            .sample(random_state=42)
            .iloc[0],
            "n_unique_5end": df["n_unique_5end"].sample(random_state=42).iloc[0],
            "n_unique_3end": df["n_unique_3end"].sample(random_state=42).iloc[0],
            "5end_gini": df["5end_gini"].sample(random_state=42).iloc[0],
            "3end_gini": df["3end_gini"].sample(random_state=42).iloc[0],
            "n_cells": df["cell_id"].nunique(),
            "n_donors": df["donor_id"].nunique(),
            "n_tissues": df["tissue"].nunique(),
            "germline_distance": df["germline_distance"].min(),
            **{l: df[l].any() for l in labels},
        }
    )


def pivot_data(feature: str, data: pd.DataFrame) -> csr_matrix:
    """
    Pivot data for a feature
    :param feature: str, feature to pivot
    :param data: pd.DataFrame, data
    """
    return (
        data.pivot_table(index="cell_id", columns="Cluster", values=feature)
        .fillna(0)
        .values
    )


class SLAVseqDataSet(object):
    def __init__(
        self,
        pqt: str,
        ref_fa: str,
        meta: pd.DataFrame,
        features: list,
        labels: list,
    ):
        self.meta = meta
        self.features = features
        self.labels = labels
        self._filtered = False

        assert Path(ref_fa).exists(), f"{ref_fa} does not exist"

        logger.info(f"Computing genomic background frequencies from {ref_fa}")
        self.ref_fa = ref_fa

        assert Path(pqt).exists(), f"{pqt} does not exist"
        self.data = pd.read_parquet(pqt)
        self.raw = self.data.copy()

        for c in [
            "Chromosome",
            "Start",
            "End",
            "locus",
            "cell_id",
            "donor_id",
            *features,
        ]:
            assert c in self.data.columns, f"Must have {c} in data"
        self.message(with_cells_donors=True)

    def message(self, with_cells_donors=False):
        n_regions, n_uniq_regions = len(self.data), self.data["locus"].nunique()
        ndonors, ncells = (
            self.data["donor_id"].nunique(),
            self.data["cell_id"].nunique(),
        )

        if with_cells_donors:
            logger.info(
                f"{n_regions} ({n_uniq_regions} unique) SLAVseq regions from {ncells} samples from {ndonors} donors"
            )
        else:
            logger.info(f"{n_regions} ({n_uniq_regions} unique) SLAVseq regions")

        if "label" in self.data.columns:
            label_counts = self.data["label"].value_counts()
            logger.info(
                f"{label_counts['KRGL']} KRGL, {label_counts['KNRGL']} KNRGL and {label_counts['OTHER']} OTHER"
            )
            avg_per_sample = (
                self.data.groupby("cell_id")["label"]
                .value_counts()
                .unstack()
                .fillna(0)
                .mean()
                .astype(int)
            )
            logger.info(
                f"Average per sample: {avg_per_sample['KRGL']} KRGL, {avg_per_sample['KNRGL']} KNRGL and {avg_per_sample['OTHER']} OTHER"
            )

    def label(self):
        logger.info("Labelling regions...")
        self.data["label"] = self.data.apply(collate_labels, axis=1)
        if not self._filtered:
            self.raw = self.data.copy()
        self.message()

    def compute_germline_distance(self, germline_df: pd.DataFrame):
        logger.info("Computing germline distance...")
        assert (
            "label" in self.data.columns
        ), "Need to label regions to compute germline distance"
        self.data = germline_distance(self.data, germline_df=germline_df)
        if not self._filtered:
            self.raw = self.data.copy()

    def cluster(self, slack: int = 0, threads: int = 1):
        """
        Cluster/merge regions, compute motif scores, and reshape to AnnData
        :param slack: int, slack to add to the end of the region
        :param threads: int, number of threads to use for cluster merging
        """

        logger.info("Clustering regions...")
        if "label" not in self.data.columns:
            logger.info("Need to label regions to merge, running label")
            self.label()
        if "Cluster" in self.data.columns:
            logger.info("Regions already clustered, overwriting...")
            self.data = self.data.drop(columns="Cluster")

        # cluster in parallel
        logger.info(f"Clustering regions with slack={slack}...")
        self.data = pr.PyRanges(self.data).cluster(slack=slack, nb_cpu=threads).df
        self.data["Cluster"] = self.data["Cluster"].astype(str)

        # merge clusters in parallel
        start = time()  # time it
        with Pool(threads) as p:
            ret_list = p.starmap_async(
                merge_clusters,
                [
                    (name, group, self.labels)
                    for name, group in self.data.groupby("Cluster")
                ],
            )
            ret_list = [r for r in ret_list.get()]
        self.locus_data = pd.DataFrame(ret_list)
        self.locus_data["Cluster"] = self.locus_data["Cluster"].astype(str)
        self.locus_data.set_index("Cluster", inplace=True)
        self.locus_data["n_tissues"] = pd.Categorical(
            self.locus_data["n_tissues"], categories=range(1, 3)
        )
        self.locus_data["label"] = self.locus_data.apply(collate_labels, axis=1)

        logger.info(
            f"Merged {self.data['locus'].nunique()} regions into {len(self.locus_data)} regions in {time() - start:.2f} seconds"
        )
        label_counts = self.locus_data["label"].value_counts()
        logger.info(
            f"{label_counts['KRGL']} KRGL, {label_counts['KNRGL']} KNRGL and {label_counts['OTHER']} OTHER"
        )

    def make_anndata(self, threads: int = 1):
        """
        Reshape data to AnnData
        :param threads: int, number of threads to use for reshaping
        """

        if getattr(self, "locus_data", None) is None:
            raise ValueError(
                "Need to cluster regions to reshape to AnnData. Please run cluster()"
            )
        if getattr(self, "ad", None) is not None:
            logger.info("AnnData already exists, skipping reshaping, overwriting...")

        start = time()
        cells, clusters = self.data["cell_id"].unique(), self.locus_data.index
        # pivot data in parallel
        with Pool(threads) as p:
            wide_data = p.starmap(pivot_data, [(f, self.data) for f in self.features])
        wide_data = {f: w for f, w in zip(self.features, wide_data)}

        self.ad = ad.AnnData(
            X=wide_data["n_reads"],
            obs=self.meta.set_index("sample_id").loc[cells],
            var=self.locus_data,
            layers={f: wide_data[f] for f in self.features if f != "n_reads"},
        )
        logger.info(f"Reshaped to AnnData in {time() - start:.2f} seconds")

    def filter(self, query):
        """
        Filter data with query, message user
        :param query: str, query to filter data
        """

        logger.info(f"Filtering data with query: {query}")
        self.data = self.data.query(query)
        self._filtered = True
        self.message()

        # filter adata too
        if getattr(self, "ad", None) is not None:
            logger.info(
                f"Filtering AnnData with {len(self.ad.obs)} cells and {len(self.ad.var)} regions..."
            )
            cells = self.ad.obs.index.isin(self.data["cell_id"].unique())
            clusters = self.ad.var.index.isin(self.data["Cluster"].unique())
            self.ad = self.ad[cells, clusters].copy()
            self.locus_data = self.locus_data.loc[self.ad.var.index]
            logger.info(
                f"Filtered AnnData to {len(self.ad.obs)} cells and {len(self.ad.var)} regions..."
            )

    def compute_motif(self, pwm, name: str, compute_bg: bool = True):
        """
        Compute motif scores for each region for a provided PWM
        :param pwm: motifs.matrix.PositionWeightMatrix, PWM to score regions
        :param name: str, name of the column to store motif scores
        :param compute_bg: bool, compute background frequencies
        """

        if getattr(self, "locus_data", None) is None:
            raise ValueError(
                "Need to cluster regions to compute motif scores. Please run cluster()"
            )

        if name in self.locus_data.columns:
            raise ValueError(f"{name} already exists in locus_data")

        if compute_bg:
            logger.info(f"Computing genomic background frequencies from {self.ref_fa}")
            bg = defaultdict(list)
            for record in SeqIO.parse(self.ref_fa, "fasta"):
                for base in "ACGT":
                    bg[base].append(record.seq.count(base) / len(record.seq))
            bg = {base: sum(b) / len(b) for base, b in bg.items()}
            pssm = pwm.log_odds(bg)
        else:
            pssm = pwm.log_odds()

        rpssm = pssm.reverse_complement()

        logger.info(f"Computing {name} motif scores for each region...")
        fa = Fasta(self.ref_fa)
        en_scores = []
        start = time()
        for d in self.locus_data.itertuples():
            chrom, start, end = d.Chromosome, d.Start, d.End
            end1 = start + 7 if end - start < 7 else end
            seq = Seq(str(fa[chrom][start:end1]))
            en_scores.append(
                max([pssm.calculate(seq).max(), rpssm.calculate(seq).max()])
            )
        logger.info(
            f"Computed {name} motif scores on {len(self.locus_data)} consensus regions in {time() - start:.2f} seconds"
        )
        self.locus_data[name] = en_scores

        if getattr(self, "ad", None) is not None:
            self.ad.var[name] = en_scores


if __name__ == "__main__":

    from scripts.get_peak_features import FEATURES

    logging.basicConfig(
        filename=snakemake.log[0],  # type: ignore
        filemode="w",
        level=logging.INFO,
    )
    logger = logging.getLogger(__name__)

    # LABELS
    LABELS = [
        "bulk",
        "megane",
        "graffite",
        "l1hs",
        "l1pa2",
        "l1pa3",
        "l1pa4",
        "l1pa5",
        "l1pa6",
        "polyA",
        "polyT",
        "primer_sites",
    ]

    # meta data
    meta = pd.read_csv(snakemake.config["samples"], sep="\t", dtype={"sample_id": str, "tissue_id": str, "donor_id": str}).drop(columns=["R1", "R2"])  # type: ignore
    donors = pd.read_csv(snakemake.config["donors"], sep="\t", dtype={"donor_id": str})  # type: ignore
    meta = meta.merge(
        donors,
        on=["donor_id", "brain_id", "libd_id", "sex", "race", "age", "diagnosis"],
    )
    meta["cell_id"] = meta["sample_id"]

    # germline
    l1hs = pr.read_bed(snakemake.input.l1hs_rmsk).df
    megane = pr.read_bed(snakemake.input.megane[0]).df
    graffite = pr.read_bed(snakemake.input.graffite[0]).df
    bulk = pd.concat(
        [pd.read_csv(f, sep="\t", low_memory=False) for f in snakemake.input.bulk]
    )
    bulk.columns = bulk.columns.str.replace("#", "")
    bulk.query("n_reads > 100 and mapq_q1 == 60", inplace=True)
    germline = pd.concat([l1hs, megane, graffite, bulk])
    germline = pr.PyRanges(germline).merge().df

    # SLAVseq object
    slavseq = SLAVseqDataSet(
        snakemake.input.pqt[0],
        ref_fa=snakemake.input.fasta,
        meta=meta,
        features=FEATURES,
        labels=LABELS,
    )
    slavseq.label()
    slavseq.compute_germline_distance(germline_df=germline)
    slavseq.cluster(threads=snakemake.threads)

    flasch = pd.read_excel(snakemake.input.flasch)
    instances = [Seq(site.replace("/", "")) for site in flasch["L1 EN Cleavage"].values]
    m = motifs.create(instances)
    pwm = m.counts.normalize()
    slavseq.compute_motif(pwm, name="flasch_score", compute_bg=True)

    cons = [Seq("TTTTAA")]
    m = motifs.create(cons)
    pwm = m.counts.normalize()
    slavseq.compute_motif(pwm, name="cons_score", compute_bg=True)

    slavseq.make_anndata(threads=snakemake.threads)

    # write object to pickle
    with open(snakemake.output[0], "wb") as f:
        pkl.dump(slavseq, f)
    print("Finished!")
