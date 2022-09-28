"""
step37.py: Correlation with sites
"""
import argparse
import itertools
import multiprocessing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import step00

input_data = pandas.DataFrame()


def pearson(a: str, b: str) -> float:
    return scipy.stats.pearsonr(input_data[a], input_data[b])[1]


def spearman(a: str, b: str) -> float:
    return scipy.stats.spearmanr(input_data[a], input_data[b])[1]


def kendall(a: str, b: str) -> float:
    return scipy.stats.kendalltau(input_data[a], input_data[b])[1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pearson", help="Pearson r", action="store_true", default=False)
    group.add_argument("--spearman", help="Spearman r", action="store_true", default=False)
    group.add_argument("--kendall", help="kendall tau", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tar.gz"):
        raise ValueError("Input file must end with .tar.gz!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUS must be a positive integer!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = step00.read_pickle(args.input)
    input_data.index = list(map(lambda x: step00.consistency_taxonomy(x, 2), list(input_data.index)))
    input_data.sort_index(inplace=True)
    input_data = input_data.groupby(input_data.index).sum().T
    print(input_data)

    sites = list(input_data.columns)

    output_data = pandas.DataFrame()
    with multiprocessing.Pool(args.cpus) as pool:
        if args.pearson:
            output_data = input_data.corr(method="pearson").stack().reset_index(name="correlation")
            output_data["-log10(p)"] = -1 * numpy.log10(pool.starmap(pearson, itertools.product(sites, repeat=2)))
        elif args.spearman:
            output_data = input_data.corr(method="spearman").stack().reset_index(name="correlation")
            output_data["-log10(p)"] = -1 * numpy.log10(pool.starmap(spearman, itertools.product(sites, repeat=2)))
        elif args.kendall:
            output_data = input_data.corr(method="kendall").stack().reset_index(name="correlation")
            output_data["-log10(p)"] = -1 * numpy.log10(pool.starmap(kendall, itertools.product(sites, repeat=2)))
        else:
            raise Exception("Something went wrong!!")
    output_data = output_data.replace(numpy.inf, 0)
    print(output_data)

    g = seaborn.relplot(data=output_data, x="level_0", y="level_1", hue="correlation", size="-log10(p)", legend="brief", height=24, aspect=1)
    g.set_xlabels("Sites")
    g.set_ylabels("Sites")

    g.savefig(args.output)
