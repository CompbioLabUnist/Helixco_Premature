"""
step84.py: Abundance/Proportion Distribution with DAT - 2
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import step00


def GW_to_float(GW):
    week, day = GW.split("+")
    return int(week) + int(day) / 7


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("DAT", type=str, help="DAT TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output TAR file")
    parser.add_argument("--site", type=str, help="Site to plot", required=True)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.DAT.endswith(".tsv"):
        raise ValueError("DAT file must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID")
    input_data["taxonomy"] = list(map(lambda x: step00.simplified_taxonomy(x) if step00.filtering_taxonomy(x) else "Unclassified", input_data["taxonomy"]))
    input_data = input_data.groupby("taxonomy").sum().T
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str, index_col="#SampleID").replace(to_replace=-1, value=None).dropna(axis="columns", how="all")
    metadata = metadata.loc[sorted(set(input_data.index) & set(metadata.index)), :]
    metadata = metadata.loc[(metadata["Site"] == args.site)]
    print(metadata)

    input_data = input_data.loc[list(filter(lambda x: x in set(input_data.index), list(metadata.index))), :]
    print(input_data)

    DAT_data = pandas.read_csv(args.DAT, sep="\t", index_col=0)
    PTB_DAT = list(map(step00.simplified_taxonomy, list(DAT_data.loc[(DAT_data["log2FoldChange"] > 1) & (DAT_data["padj"] < 0.05)].sort_values("log2FoldChange").index)))
    Normal_DAT = list(map(step00.simplified_taxonomy, list(DAT_data.loc[(DAT_data["log2FoldChange"] < -1) & (DAT_data["padj"] < 0.05)].sort_values("log2FoldChange").index)))
    print(DAT_data)

    taxa = Normal_DAT + PTB_DAT + sorted(list(filter(lambda x: (x not in Normal_DAT) and (x not in PTB_DAT), list(input_data.columns))), key=lambda x: sum(input_data[x]), reverse=True)
    input_data = input_data.loc[sorted(list(input_data.index), key=lambda x: (sum(input_data.loc[x, Normal_DAT]) - sum(input_data.loc[x, PTB_DAT]), metadata.loc[x, "Gestational Week"]), reverse=True), taxa]
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    output_data = pandas.DataFrame(index=input_data.index)
    output_data["Normal_DAT"] = numpy.sum(input_data.loc[:, Normal_DAT], axis=1)
    output_data["PTB_DAT"] = numpy.sum(input_data.loc[:, PTB_DAT], axis=1)
    output_data["PTB-Normal"] = output_data["PTB_DAT"] - output_data["Normal_DAT"]
    output_data["GW"] = list(map(lambda x: GW_to_float(metadata.loc[x, "Detail Gestational Week"]), list(output_data.index)))
    print(output_data)

    r, p = scipy.stats.pearsonr(output_data["GW"], output_data["PTB-Normal"])
    print(r, p)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
    seaborn.regplot(data=output_data, x="GW", y="PTB-Normal", scatter_kws={"linewidth": 0, "s": 400}, ax=ax)
    matplotlib.pyplot.axvline(x=37, color="k", linestyle="--")

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(args.output)
    matplotlib.pyplot.close(fig)
