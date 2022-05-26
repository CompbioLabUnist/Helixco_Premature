"""
step50.py: Volcano plot with differetially expressed taxa
"""
import argparse
import multiprocessing
import tarfile
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


def ratio(taxon: str, site: str) -> float:
    try:
        return numpy.log2(numpy.mean(input_data.loc[(input_data["Simple Premature"] == "Early PTB") & (input_data["Site"] == site), pandas.MultiIndex.from_tuples([taxon])].to_numpy()) / numpy.mean(input_data.loc[(input_data["Simple Premature"] == "Late PTB+Normal") & (input_data["Site"] == site), pandas.MultiIndex.from_tuples([taxon])].to_numpy()))
    except RuntimeWarning:
        return 0.0


def pvalue(taxon: str, site: str) -> float:
    try:
        return -1 * numpy.log10(scipy.stats.mannwhitneyu(input_data.loc[(input_data["Simple Premature"] == "Early PTB") & (input_data["Site"] == site), pandas.MultiIndex.from_tuples([taxon])].to_numpy(), input_data.loc[(input_data["Simple Premature"] == "Late PTB+Normal") & (input_data["Site"] == site), pandas.MultiIndex.from_tuples([taxon])].to_numpy(), alternative="two-sided")[1])
    except ValueError:
        return 1.0


def multiply(ratio: float, pvalue: float) -> float:
    return abs(ratio * pvalue)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tsv file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    if not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata must end with .tsv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be positive!!")

    input_data = pandas.read_csv(args.input, skiprows=1, sep="\t", index_col=["taxonomy", "#OTU ID"]).T
    for index in tqdm.tqdm(list(input_data.index)):
        input_data.loc[index, :] = input_data.loc[index, :] / sum(input_data.loc[index, :])
    taxa = input_data.columns
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col=0)
    print(metadata)
    print(sorted(metadata.columns))

    input_data = pandas.concat([input_data, metadata], axis="columns", verify_integrity=True, join="inner")
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    warnings.filterwarnings("error")

    figures = list()

    for site in tqdm.tqdm(step00.selected_long_sites):
        figures.append(f"{site}.pdf")
        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
        texts = list()

        output_data = pandas.DataFrame(index=taxa, dtype=float)
        print(output_data)

        with multiprocessing.Pool(args.cpus) as pool:
            output_data["log2(EP/LP+F)"] = pool.starmap(ratio, [(taxon, site) for taxon in taxa])
            output_data["-log10(p)"] = pool.starmap(pvalue, [(taxon, site) for taxon in taxa])
            output_data["Multiple"] = pool.starmap(multiply, zip(output_data["log2(EP/LP+F)"], output_data["-log10(p)"]))
        output_data.sort_values("Multiple", ascending=False, inplace=True)
        print(output_data)

        down_results = output_data.loc[((output_data["log2(EP/LP+F)"] < numpy.log2(1 / ratio_threshold)) & (output_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
        up_results = output_data.loc[((output_data["log2(EP/LP+F)"] > numpy.log2(ratio_threshold)) & (output_data["-log10(p)"] > (-1 * numpy.log10(p_threshold)))), :]
        ns_results = output_data.loc[(((output_data["log2(EP/LP+F)"] < numpy.log2(ratio_threshold)) & (output_data["log2(EP/LP+F)"] > numpy.log2(1 / ratio_threshold))) | (output_data["-log10(p)"] < (-1 * numpy.log10(p_threshold)))), :]

        matplotlib.pyplot.scatter(ns_results["log2(EP/LP+F)"], ns_results["-log10(p)"], s=100, c="gray", marker="o", edgecolors=None, label="NS")
        matplotlib.pyplot.scatter(up_results["log2(EP/LP+F)"], up_results["-log10(p)"], s=100, c="red", marker="o", edgecolors=None, label="Up")
        matplotlib.pyplot.scatter(down_results["log2(EP/LP+F)"], down_results["-log10(p)"], s=100, c="blue", marker="o", edgecolors=None, label="Down")

        for index, row in down_results.iloc[:5, :].iterrows():
            texts.append(matplotlib.pyplot.text(row["log2(EP/LP+F)"], row["-log10(p)"], step00.consistency_taxonomy(index[0], 2), color="black", fontsize="xx-small"))

        for index, row in up_results.iloc[:5, :].iterrows():
            texts.append(matplotlib.pyplot.text(row["log2(EP/LP+F)"], row["-log10(p)"], step00.consistency_taxonomy(index[0], 2), color="black", fontsize="xx-small"))

        matplotlib.pyplot.xlabel("log2(EP/LP+F)")
        matplotlib.pyplot.ylabel("-log10(p)")
        matplotlib.pyplot.title(site)
        matplotlib.pyplot.axvline(numpy.log2(1 / ratio_threshold), color="k", linestyle="--")
        matplotlib.pyplot.axvline(numpy.log2(ratio_threshold), color="k", linestyle="--")
        matplotlib.pyplot.axhline(-1 * numpy.log10(p_threshold), color="k", linestyle="--")
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.tight_layout()

        adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black", alpha=0.3), lim=step00.small, ax=ax)
        matplotlib.pyplot.tight_layout()

        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
