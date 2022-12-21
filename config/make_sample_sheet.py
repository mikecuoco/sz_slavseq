#!/usr/bin/env python
# Created on: 10/26/22, 1:59 PM
__author__ = "Michael Cuoco"

import pandas as pd
import gzip, re, subprocess, argparse

# define funtions for file processing
def my_files(dirname):
    """yield fastq.gz files in a directory"""
    command = 'find {} -name "*.fastq.gz"'.format(dirname)
    with subprocess.Popen(command, shell=True, stdout=subprocess.PIPE) as p:
        for line in p.stdout:
            l = line.decode().rstrip()
            yield l


def num_reads(file):
    """get number of reads in an open file"""
    return int((1 + sum(1 for _ in file)) / 4)


def file_info(filename):
    """get info for a gzipped file"""
    first = None
    with gzip.open(filename, "rt") as f:
        first = f.readline().rstrip()
        nreads = num_reads(f)

    return (filename, first, nreads)


# define functions for file meta data generate from file name
def fields(filename):
    """parse filename to get relevant fields"""
    m = re.search("/([^/]+?)(_001)?.fastq.gz$", filename)
    assert m is not None
    basename = m.group(1)

    r = dict()
    r["filename"] = filename
    r["basename"] = basename
    r["pair_id"] = re.sub("_(R[12])$", "", basename)

    if basename.startswith("plate"):
        m2 = re.search("^(plate\d+)_([A-H]\d+)_(S\d+)_(R[12])$", basename)
        assert m2 is not None
        r["individual"] = "CommonBrain"
        r["sample_id1"] = m2.group(1).lower() + m2.group(2).upper()
        r["sample_id2"] = m2.group(3).upper()
        r["read"] = m2.group(4).upper()
        r["tissue"] = "brain"
        r["tissue_id"] = "CommonBrain"
        r["sample_type"] = "single_cell"
        r["dna_type"] = "mda"
    elif basename.upper().startswith("US"):
        m2 = re.search(
            "^(US([DH])(\d+))_?([A-H]\d+)_(S\d+)_(R[12])$",
            basename,
            flags=re.IGNORECASE,
        )
        r["individual"] = re.sub("^0+", "", m2.group(3))
        r["sample_id1"] = m2.group(4).upper()
        r["sample_id2"] = m2.group(5).upper()
        r["read"] = m2.group(6).upper()
        r["tissue"] = "HIPPO" if m2.group(2).upper() == "H" else "DLPFC"
        ind = r["individual"] if int(r["individual"]) >= 10 else "0" + r["individual"]
        r["tissue_id"] = "US" + m2.group(2).upper() + ind
        r["sample_type"] = "single_cell"
        r["dna_type"] = "mda"
    elif basename.startswith("gDNA"):
        m2 = re.search("^gDNA_(US([DH])(\d+))_(R[12])$", basename, flags=re.IGNORECASE)
        r["individual"] = re.sub("^0+", "", m2.group(3))
        r["sample_id1"] = "bulk"
        r["sample_id2"] = "Sbulk"
        r["read"] = m2.group(4).upper()
        r["tissue"] = "HIPPO" if m2.group(2).upper() == "H" else "DLPFC"
        ind = r["individual"] if int(r["individual"]) >= 10 else "0" + r["individual"]
        r["tissue_id"] = "US" + m2.group(2).upper() + ind
        r["sample_type"] = "bulk"
        r["dna_type"] = "bulk"
    else:
        r = None

    return r


def print_unique(df):
    print(
        f"""
		{df.shape[0]} files
		{len(set(df['first_read']))} unique files based on first read
		{len(set(df['first_read'] + str(df['num_reads'])))} unique files baesd on first read + num reads
		"""
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dirs", nargs="+", help="directories to search for fastq.gz files")
    parser.add_argument("--meta", type=str, help="path to slavseq_metadata.csv")
    parser.add_argument("--find-unique", action="store_true", help="find unique files", default=False)
    parser.add_argument("--threads", type=int, nargs=1, help="# threads to use for unique file finding", default=1)
    args = parser.parse_args()

    if args.find_unique:
        from joblib import Parallel, delayed
        print("searching for fastq.gz files in {}".format(", ".join(args.data_dirs)))
        # use threads because first_line is i/o bound
        results = Parallel(n_jobs=args.threads)(
            delayed(file_info)(f) for d in args.data_dirs for f in my_files(d)
        )
        print("found {} files".format(len(results)))

        # convert to dataframe
        print("extracting first read from each file")
        logg_files = pd.DataFrame.from_records(
            results, columns=["filename", "first_read", "num_reads"]
        )

        # remove duplicates
        print("removing duplicate fastqs")
        logg_files.drop_duplicates(
            subset=["first_read", "num_reads"], inplace=True
        ) 
    else:
        logg_files = pd.DataFrame({"filename": [f for d in args.data_dirs for f in my_files(d)]})
    
    # extract sample metadata from filenames
    df = pd.DataFrame.from_records([dict(fields(x)) for x in logg_files["filename"]]).merge(logg_files, on="filename")

    meta = pd.read_csv(args.meta, sep="\t").drop(["MDA_PERFORMED","BULK_PERFORMED"], axis = 1)
    meta.columns = meta.columns.str.lower()
    meta.set_index("tissue_id", inplace=True)

    # write to tsv
    (
        df.filter(
            items=[
                "filename",
                "read",
                "individual",
                "pair_id",
                "dna_type",
                "tissue_id",
            ],
            axis=1,
        )  # keep only relevant columns
        .pivot(
            columns="read",
            values="filename",
            index=["pair_id", "individual", "dna_type", "tissue_id"],
        )  # pivot to get R1 and R2 in same row
        .reset_index()
        .join(meta, on="tissue_id", how="left")  # join with metadata
        .rename(columns={"pair_id": "sample", "individual": "donor"})
        .sort_values("donor")  # sort by donor_id
        .to_csv("all_donors.tsv", sep="\t", index=False)
    )

df = pd.read_csv('all_donors.tsv', sep='\t')
df[df["donor"].isin(["2","20","37","28"])].to_csv("four_donors.tsv", sep='\t', index=False)