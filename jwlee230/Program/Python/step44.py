"""
step44.py: Volcano plot with differetially expressed taxa
"""
import argparse
import warnings
import adjustText
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import step00

ratio_threshold = 2
p_threshold = 0.05


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    input_data["-log10(p)"] = -1 * numpy.log10(input_data["padj"])
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    warnings.filterwarnings("error")

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    texts = list()

    down_results = input_data.loc[((input_data["log2FoldChange"] < numpy.log2(1 / ratio_threshold)) & (input_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
    up_results = input_data.loc[((input_data["log2FoldChange"] > numpy.log2(ratio_threshold)) & (input_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
    ns_results = input_data.loc[(((input_data["log2FoldChange"] < numpy.log2(ratio_threshold)) & (input_data["log2FoldChange"] > numpy.log2(1 / ratio_threshold))) | (input_data["-log10(p)"] < (-1 * numpy.log10(p_threshold)))), :]

    matplotlib.pyplot.scatter(ns_results["log2FoldChange"], ns_results["-log10(p)"], s=100, c="gray", marker="o", edgecolors=None)
    matplotlib.pyplot.scatter(up_results["log2FoldChange"], up_results["-log10(p)"], s=100, c="red", marker="o", edgecolors=None)
    matplotlib.pyplot.scatter(down_results["log2FoldChange"], down_results["-log10(p)"], s=100, c="blue", marker="o", edgecolors=None)

    for index, row in down_results.iterrows():
        print("Down /", index)
        texts.append(matplotlib.pyplot.text(row["log2FoldChange"], row["-log10(p)"], step00.simplified_taxonomy(index), color="black", fontsize="small"))

    for index, row in up_results.iterrows():
        print("Up /", index)
        texts.append(matplotlib.pyplot.text(row["log2FoldChange"], row["-log10(p)"], step00.simplified_taxonomy(index), color="black", fontsize="small"))

    matplotlib.pyplot.xlabel("log2FoldChange")
    matplotlib.pyplot.ylabel("-log10(p)")
    matplotlib.pyplot.axvline(numpy.log2(1 / ratio_threshold), color="k", linestyle="--")
    matplotlib.pyplot.axvline(numpy.log2(ratio_threshold), color="k", linestyle="--")
    matplotlib.pyplot.axhline(-1 * numpy.log10(p_threshold), color="k", linestyle="--")
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.tight_layout()

    adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black", alpha=0.3), lim=step00.small, ax=ax)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
