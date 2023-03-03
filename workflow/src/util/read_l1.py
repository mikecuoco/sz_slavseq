import pandas as pd
import pyranges as pr


def read_rmsk(fn):
    # read the rmsk file
    rmsk = pd.read_csv(
        fn,
        skiprows=3,
        delim_whitespace=True,
        names=["Chromosome", "Start", "End", "Strand", "repeat"],
        usecols=[4, 5, 6, 8, 9],
    )

    # filter for rep_names
    rep_names = [
        "L1HS_3end",
        "L1PA2_3end",
        "L1PA3_3end",
        "L1PA4_3end",
        "L1PA5_3end",
        "L1PA6_3end",
    ]
    assert all(
        [x in rmsk["repeat"].unique() for x in rep_names]
    ), "Not all rep_names are in rmsk"
    rmsk = rmsk[rmsk["repeat"].isin(rep_names)]

    # correct strand
    rmsk["Strand"] = rmsk.apply(lambda x: "+" if x.Strand == "+" else "-", axis=1)

    return rmsk


def read_knrgl(fn):

    # read the knrgl bed file
    knrgl = pd.read_csv(
        fn,
        sep="\t",
        header=None,
        names=["Chromosome", "Start", "End", "Strand", "SVLEN", "SVTYPE"],
        dtype={"Chromosome": str, "Start": int, "End": int},
    )
    assert knrgl.shape[1] == 6, "knrgl bed file has wrong number of columns"

    # adjust start and end
    knrgl["Start"] = knrgl.apply(
        lambda x: x.Start - 750 if x.Strand == "-" else x.Start, axis=1
    )
    knrgl["End"] = knrgl.apply(
        lambda x: x.End + 750 if x.Strand == "+" else x.End, axis=1
    )

    return knrgl
