"""
step83.py: Abundance/Proportion Distribution with DAT
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
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
    input_data = input_data.loc[sorted(list(input_data.index), key=lambda x: (sum(input_data.loc[x, PTB_DAT]) - sum(input_data.loc[x, Normal_DAT]), metadata.loc[x, "Gestational Week"]), reverse=True), taxa]
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, axs = matplotlib.pyplot.subplots(figsize=(32, 18), nrows=2, gridspec_kw={"height_ratios": [2, 1]})

    axs[0].bar(range(input_data.shape[0]), numpy.sum(input_data.loc[:, Normal_DAT], axis=1) * -1, color=step00.PTB_two_colors["Normal"], align="edge", linewidth=0, label="Normal-enriched DAT")
    axs[0].bar(range(input_data.shape[0]), numpy.sum(input_data.loc[:, PTB_DAT], axis=1), color=step00.PTB_two_colors["PTB"], align="edge", linewidth=0, label="PTB-enriched DAT")
    axs[0].plot(numpy.array(range(input_data.shape[0])) + 0.5, numpy.sum(input_data.loc[:, PTB_DAT], axis=1) - numpy.sum(input_data.loc[:, Normal_DAT], axis=1), color="black", linestyle="--", linewidth=4, label="PTB DAT - Normal DAT")

    axs[0].set_xticks([])
    axs[0].set_xlabel("")
    axs[0].set_ylabel("DAT proportion")
    axs[0].legend(loc="upper center", fontsize="xx-small")
    axs[0].grid(True)

    axs[1].bar(range(input_data.shape[0]), list(map(lambda x: GW_to_float(metadata.loc[x, "Detail Gestational Week"]) if (metadata.loc[x, "Premature"] == "PTB") else 0, list(input_data.index))), color="tab:red", linewidth=0, label="PTB subject")
    axs[1].bar(range(input_data.shape[0]), list(map(lambda x: GW_to_float(metadata.loc[x, "Detail Gestational Week"]) if (metadata.loc[x, "Premature"] == "Normal") else 0, list(input_data.index))), color="tab:green", linewidth=0, label="Normal subject")
    axs[1].scatter([i for i in range(input_data.shape[0]) if (metadata.loc[list(input_data.index)[i], "PROM"] == "True")], [40 for _ in list(filter(lambda x: x == "True", metadata["PROM"]))], marker="^", c="black", edgecolors=None, label="PROM")

    axs[1].set_xticks([])
    axs[1].set_xlabel(f"{input_data.shape[0]} samples")
    axs[1].set_ylabel("GW")
    axs[1].set_ylim(20, 45)
    axs[1].legend(loc="lower left")
    axs[1].grid(True)
    axs[1].axhline(37, linestyle="--", color="black", linewidth=4)

    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
