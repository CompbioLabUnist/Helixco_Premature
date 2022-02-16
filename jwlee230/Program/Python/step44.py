"""
step44.py: Volcano plot with differetially expressed taxa
"""
import argparse
import multiprocessing
import warnings
import adjustText
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import tqdm
import step00

input_data = pandas.DataFrame()
ratio_threshold = 2
p_threshold = 0.05


def ratio(taxo: str) -> float:
    try:
        return numpy.log2(numpy.mean(input_data.loc[(input_data["Detail Premature"] == "Early PTB"), taxo]) / numpy.mean(input_data.loc[(input_data["Detail Premature"] == "Normal"), taxo]))
    except RuntimeWarning:
        return 0


def pvalue(taxo: str) -> float:
    try:
        return -1 * numpy.log10(scipy.stats.mannwhitneyu(input_data.loc[(input_data["Detail Premature"] == "Early PTB"), taxo], input_data.loc[(input_data["Detail Premature"] == "Normal"), taxo], alternative="two-sided")[1])
    except ValueError:
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)

    args = parser.parse_args()

    if not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata must end with .tsv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    input_data = step00.read_pickle(args.input)
    input_data.index = list(map(step00.simplified_taxonomy, list(input_data.index)))
    input_data.sort_index(inplace=True)
    input_data = input_data.groupby(input_data.index).sum().T
    taxa = list(input_data.columns)
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col=0)
    print(metadata)
    print(sorted(metadata.columns))

    input_data = pandas.concat([input_data, metadata], axis="columns", verify_integrity=True, join="inner")
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    warnings.filterwarnings("error")

    output_data = pandas.DataFrame(data=numpy.zeros((len(taxa), 2)), index=taxa, columns=["log2(EP/F)", "-log10(p)"])

    with multiprocessing.Pool(args.cpus) as pool:
        output_data.loc[:, "log2(EP/F)"] = pool.map(ratio, taxa)
        output_data.loc[:, "-log10(p)"] = pool.map(pvalue, taxa)
    print(output_data)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    texts = list()

    down_results = output_data.loc[((output_data["log2(EP/F)"] < numpy.log2(1 / ratio_threshold)) & (output_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
    up_results = output_data.loc[((output_data["log2(EP/F)"] > numpy.log2(ratio_threshold)) & (output_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
    ns_results = output_data.loc[(((output_data["log2(EP/F)"] < numpy.log2(ratio_threshold)) & (output_data["log2(EP/F)"] > numpy.log2(1 / ratio_threshold))) | (output_data["-log10(p)"] < (-1 * numpy.log10(p_threshold)))), :]

    matplotlib.pyplot.scatter(ns_results["log2(EP/F)"], ns_results["-log10(p)"], s=100, c="gray", marker="o", edgecolors=None, label="NS")
    matplotlib.pyplot.scatter(up_results["log2(EP/F)"], up_results["-log10(p)"], s=100, c="red", marker="o", edgecolors=None, label="Up")
    matplotlib.pyplot.scatter(down_results["log2(EP/F)"], down_results["-log10(p)"], s=100, c="blue", marker="o", edgecolors=None, label="Down")

    for index, row in tqdm.tqdm(down_results.iterrows()):
        texts.append(matplotlib.pyplot.text(row["log2(EP/F)"], row["-log10(p)"], step00.consistency_taxonomy(index, 1), color="black", fontsize="xx-small"))

    for index, row in tqdm.tqdm(up_results.iterrows()):
        texts.append(matplotlib.pyplot.text(row["log2(EP/F)"], row["-log10(p)"], step00.consistency_taxonomy(index, 1), color="black", fontsize="xx-small"))

    matplotlib.pyplot.xlabel("log2(EP/F)")
    matplotlib.pyplot.ylabel("-log10(p)")
    matplotlib.pyplot.axvline(numpy.log2(1 / ratio_threshold), color="k", linestyle="--")
    matplotlib.pyplot.axvline(numpy.log2(ratio_threshold), color="k", linestyle="--")
    matplotlib.pyplot.axhline(-1 * numpy.log10(p_threshold), color="k", linestyle="--")
    matplotlib.pyplot.legend()
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.tight_layout()

    adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black", alpha=0.3), lim=step00.small, ax=ax)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
