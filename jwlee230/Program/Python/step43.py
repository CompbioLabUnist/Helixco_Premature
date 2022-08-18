"""
step43.py: Volcano plot with Corrleation
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
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
correlation_threshold = 0.3
p_threshold = 0.05


def pearson(clinical: str, taxonomy: str, site: str) -> typing.Tuple[float, float]:
    drawing_data = input_data.loc[((input_data[clinical] != -1) & (input_data["Site"] == site)), [taxonomy, clinical]]
    try:
        return scipy.stats.pearsonr(drawing_data[taxonomy], drawing_data[clinical])
    except scipy.stats.PearsonRConstantInputWarning:
        return (0, 1)
    except scipy.stats.PearsonRNearConstantInputWarning:
        return (0, 1)


def spearman(clinical: str, taxonomy: str, site: str) -> typing.Tuple[float, float]:
    drawing_data = input_data.loc[((input_data[clinical] != -1) & (input_data["Site"] == site)), [taxonomy, clinical]]
    return scipy.stats.spearmanr(drawing_data[taxonomy], drawing_data[clinical])


def kendall(clinical: str, taxonomy: str, site: str) -> typing.Tuple[float, float]:
    drawing_data = input_data.loc[((input_data[clinical] != -1) & (input_data["Site"] == site)), [taxonomy, clinical]]
    return scipy.stats.kendalltau(drawing_data[taxonomy], drawing_data[clinical])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pearson", help="Pearson r", action="store_true", default=False)
    group.add_argument("--spearman", help="Spearman r", action="store_true", default=False)
    group.add_argument("--kendall", help="kendall tau", action="store_true", default=False)

    args = parser.parse_args()

    if not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata must end with .tsv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")
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

    figures = list()
    with multiprocessing.Pool(args.cpus) as pool:
        for clinical, site in tqdm.tqdm(list(itertools.product(step00.numeric_columns, step00.selected_long_sites))):
            figures.append("{1}_{0}.pdf".format(clinical.replace(" ", "_"), site))
            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
            texts = list()

            if args.pearson:
                results = list(map(lambda x: (x[0][0], x[0][1], x[1]), zip(pool.starmap(pearson, [(clinical, taxo, site) for taxo in taxa]), taxa)))
            elif args.spearman:
                results = list(map(lambda x: (x[0][0], x[0][1], x[1]), zip(pool.starmap(pearson, [(clinical, taxo, site) for taxo in taxa]), taxa)))
            elif args.kendall:
                results = list(map(lambda x: (x[0][0], x[0][1], x[1]), zip(pool.starmap(pearson, [(clinical, taxo, site) for taxo in taxa]), taxa)))
            else:
                raise Exception("Something went wrong!!")

            results = list(map(lambda x: (x[0], -1 * numpy.log10(x[1]), x[2]), results))

            down_results = list(filter(lambda x: (x[0] < (-1 * correlation_threshold)) and (x[1] > -1 * numpy.log10(p_threshold)), results))
            up_results = list(filter(lambda x: (x[0] > correlation_threshold) and (x[1] > -1 * numpy.log10(p_threshold)), results))
            ns_results = list(filter(lambda x: ((-1 * correlation_threshold) < x[0] < correlation_threshold) or (x[1] < -1 * numpy.log10(p_threshold)), results))

            matplotlib.pyplot.scatter(list(map(lambda x: x[0], ns_results)), list(map(lambda x: x[1], ns_results)), s=100, c="gray", marker="o", edgecolors=None, label="NS")
            matplotlib.pyplot.scatter(list(map(lambda x: x[0], down_results)), list(map(lambda x: x[1], down_results)), s=100, c="blue", marker="o", edgecolors=None, label="Negative")
            matplotlib.pyplot.scatter(list(map(lambda x: x[0], up_results)), list(map(lambda x: x[1], up_results)), s=100, c="red", marker="o", edgecolors=None, label="Positive")

            for x, y, text in down_results + up_results:
                texts.append(matplotlib.pyplot.text(x, y, step00.consistency_taxonomy(text, 1), color="black", fontsize="xx-small"))

            matplotlib.pyplot.xlabel("Correlation")
            matplotlib.pyplot.ylabel("-log10(p)")
            matplotlib.pyplot.title("{0} in {1}".format(clinical, site))
            matplotlib.pyplot.xlim(-1, 1)
            matplotlib.pyplot.axvline(-1 * correlation_threshold, color="k", linestyle="--")
            matplotlib.pyplot.axvline(correlation_threshold, color="k", linestyle="--")
            matplotlib.pyplot.axhline(-1 * numpy.log10(p_threshold), color="k", linestyle="--")
            matplotlib.pyplot.legend()
            matplotlib.pyplot.grid(True)
            matplotlib.pyplot.tight_layout()

            adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black", alpha=0.3), lim=step00.small, ax=ax)

            fig.savefig(figures[-1])
            matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for file in tqdm.tqdm(figures):
            tar.add(file, arcname=file)
