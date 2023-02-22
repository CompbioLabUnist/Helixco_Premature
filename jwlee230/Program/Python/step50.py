"""
step50.py: Volcano plot with differetially expressed taxa
"""
import argparse
import multiprocessing
import adjustText
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import step00

input_data = pandas.DataFrame()
ratio_threshold = 2
p_threshold = 0.05


def multiply(ratio: float, pvalue: float) -> float:
    return abs(ratio * pvalue)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tsv file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output pdf file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    if not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata must end with .tsv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with pdf!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    input_data = pandas.read_csv(args.input, sep="\t", names=["taxonomy", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"], index_col="taxonomy", skiprows=1)
    input_data["-log10(p)"] = -1 * numpy.log10(input_data["padj"])
    input_data.dropna(axis="index", inplace=True)
    with multiprocessing.Pool(args.cpus) as pool:
        input_data["multiple"] = pool.starmap(multiply, zip(input_data["log2FoldChange"], input_data["-log10(p)"]))
    input_data.sort_values("multiple", ascending=False, inplace=True)
    taxa = list(filter(step00.filtering_taxonomy, list(input_data.index)))
    input_data = input_data.loc[taxa, :]
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col=0)
    print(metadata)
    print(sorted(metadata.columns))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    texts = list()

    down_results = input_data.loc[((input_data["log2FoldChange"] < numpy.log2(1 / ratio_threshold)) & (input_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
    up_results = input_data.loc[((input_data["log2FoldChange"] > numpy.log2(ratio_threshold)) & (input_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
    ns_results = input_data.loc[(((input_data["log2FoldChange"] < numpy.log2(ratio_threshold)) & (input_data["log2FoldChange"] > numpy.log2(1 / ratio_threshold))) | (input_data["-log10(p)"] < (-1 * numpy.log10(p_threshold)))), :]

    matplotlib.pyplot.scatter(ns_results["log2FoldChange"], ns_results["-log10(p)"], s=100, c="gray", marker="o", edgecolors=None, label="NS")
    matplotlib.pyplot.scatter(up_results["log2FoldChange"], up_results["-log10(p)"], s=100, c="red", marker="o", edgecolors=None, label="Up")
    matplotlib.pyplot.scatter(down_results["log2FoldChange"], down_results["-log10(p)"], s=100, c="blue", marker="o", edgecolors=None, label="Down")

    for index, row in down_results.iterrows():
        texts.append(matplotlib.pyplot.text(row["log2FoldChange"], row["-log10(p)"], step00.simplified_taxonomy(index), color="black", fontsize="xx-small"))

    for index, row in up_results.iterrows():
        texts.append(matplotlib.pyplot.text(row["log2FoldChange"], row["-log10(p)"], step00.simplified_taxonomy(index), color="black", fontsize="xx-small"))

    matplotlib.pyplot.xlabel("log2(EP/LP+F)")
    matplotlib.pyplot.ylabel("-log10(p)")
    matplotlib.pyplot.axvline(numpy.log2(1 / ratio_threshold), color="k", linestyle="--")
    matplotlib.pyplot.axvline(numpy.log2(ratio_threshold), color="k", linestyle="--")
    matplotlib.pyplot.axhline(-1 * numpy.log10(p_threshold), color="k", linestyle="--")
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.tight_layout()

    adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black", alpha=0.3), lim=step00.small, ax=ax)
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
