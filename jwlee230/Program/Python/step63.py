"""
step63.py: Volcano plot with differetially abundant taxa
"""
import argparse
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
    input_data["-log10(p)"] = -1 * numpy.log10(input_data["pvalues"])
    print(input_data)

    input_data = input_data.loc[list(map(step00.filtering_taxonomy, list(input_data.index))), ].sort_values("-log10(p)", ascending=False)
    input_data["simple_name"] = list(map(step00.consistency_taxonomy, list(input_data.index)))
    print(input_data)

    ceil = numpy.ceil(max(numpy.absolute(input_data["logFC"])))
    input_data.to_csv(args.output.replace(".pdf", ".list.tsv"), sep="\t")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    texts = list()

    down_results = input_data.loc[((input_data["logFC"] < numpy.log2(1 / ratio_threshold)) & (input_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
    up_results = input_data.loc[((input_data["logFC"] > numpy.log2(ratio_threshold)) & (input_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
    ns_results = input_data.loc[(((input_data["logFC"] < numpy.log2(ratio_threshold)) & (input_data["logFC"] > numpy.log2(1 / ratio_threshold))) | (input_data["-log10(p)"] < (-1 * numpy.log10(p_threshold)))), :]

    print(up_results.sort_values("simple_name"))
    print(down_results.sort_values("simple_name"))

    matplotlib.pyplot.scatter(ns_results["logFC"], ns_results["-log10(p)"], s=100, c="gray", marker="o", edgecolors=None)
    matplotlib.pyplot.scatter(up_results["logFC"], up_results["-log10(p)"], s=100, c="red", marker="o", edgecolors=None)
    matplotlib.pyplot.scatter(down_results["logFC"], down_results["-log10(p)"], s=100, c="blue", marker="o", edgecolors=None)

    matplotlib.pyplot.text(numpy.log2(1 / ratio_threshold), 0, f"log2(FC)={numpy.log2(1 / ratio_threshold)}", c="k", fontsize="xx-small", horizontalalignment="right", verticalalignment="baseline", rotation="vertical")
    matplotlib.pyplot.text(numpy.log2(ratio_threshold), 0, f"log2(FC)={numpy.log2(ratio_threshold)}", c="k", fontsize="xx-small", horizontalalignment="right", verticalalignment="baseline", rotation="vertical")
    matplotlib.pyplot.text(-ceil, -1 * numpy.log10(p_threshold), f"p={p_threshold:.2f}", c="k", fontsize="xx-small", horizontalalignment="left", verticalalignment="baseline")

    for index, row in down_results.iloc[:5, :].iterrows():
        texts.append(matplotlib.pyplot.text(row["logFC"], row["-log10(p)"], f"{step00.simplified_taxonomy(index)}", color="black", fontsize="small"))

    for index, row in up_results.iloc[:5, :].iterrows():
        texts.append(matplotlib.pyplot.text(row["logFC"], row["-log10(p)"], f"{step00.simplified_taxonomy(index)}", color="black", fontsize="small"))

    matplotlib.pyplot.xlabel("logFC")
    matplotlib.pyplot.ylabel("-log10(p)")
    matplotlib.pyplot.axvline(numpy.log2(1 / ratio_threshold), color="k", linestyle="--")
    matplotlib.pyplot.axvline(numpy.log2(ratio_threshold), color="k", linestyle="--")
    matplotlib.pyplot.axhline(-1 * numpy.log10(p_threshold), color="k", linestyle="--")
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.title(f"Up: {up_results.shape[0]}; Down: {down_results.shape[0]}")
    matplotlib.pyplot.xlim((-ceil, ceil))
    matplotlib.pyplot.xlabel("log2(FoldChange)")
    matplotlib.pyplot.tight_layout()

    adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black", alpha=0.3), lim=step00.big, ax=ax)
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
