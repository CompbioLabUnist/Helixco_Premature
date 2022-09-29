"""
step36-1.py: Draw cluster map plot within bacteria for Macrogen data
"""
import argparse
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
    parser.add_argument("output", help="Output PNG file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)
    parser.add_argument("--correlation", help="Correlation threshold", type=float, default=0.8)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.001)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pearson", help="Pearson r", action="store_true", default=False)
    group.add_argument("--spearman", help="Spearman r", action="store_true", default=False)
    group.add_argument("--kendall", help="kendall tau", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tar.gz"):
        raise ValueError("Input file must end with .tar.gz!!")
    elif not args.output.endswith(".png"):
        raise ValueError("Output file must end with .PNG!!")
    elif args.cpus < 1:
        raise ValueError("CPUS must be a positive integer!!")
    elif not (0 < args.correlation < 1):
        raise ValueError("Correlation must be (0, 1)")
    elif not (0 < args.p < 1):
        raise ValueError("P-value must be (0, 1)")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = step00.read_pickle(args.input).T
    input_data.columns = list(map(step00.simplified_taxonomy, list(input_data.columns)))
    input_data = input_data.iloc[1:, :]
    print(input_data)

    pathogens = list(input_data.columns)
    print("Pathgens:", len(pathogens))

    output_data = pandas.DataFrame()
    with multiprocessing.Pool(args.cpus) as pool:
        if args.pearson:
            output_data = input_data.corr(method="pearson").stack().reset_index(name="correlation")
            output_data["-log10(p)"] = -1 * numpy.log10(pool.starmap(pearson, output_data[["level_0", "level_1"]].itertuples(index=False, name=None)))
        elif args.spearman:
            output_data = input_data.corr(method="spearman").stack().reset_index(name="correlation")
            output_data["-log10(p)"] = -1 * numpy.log10(pool.starmap(spearman, output_data[["level_0", "level_1"]].itertuples(index=False, name=None)))
        elif args.kendall:
            output_data = input_data.corr(method="kendall").stack().reset_index(name="correlation")
            output_data["-log10(p)"] = -1 * numpy.log10(pool.starmap(kendall, output_data[["level_0", "level_1"]].itertuples(index=False, name=None)))
        else:
            raise Exception("Something went wrong!!")
    output_data = output_data.replace(numpy.inf, numpy.nan).dropna(subset=["-log10(p)"], axis="index")
    output_data = output_data.loc[(numpy.abs(output_data["correlation"]) > args.correlation) & (output_data["-log10(p)"] > -1 * numpy.log10(args.p))]
    output_data["x"] = list(map(lambda x: pathogens.index(x), output_data["level_0"]))
    output_data["y"] = list(map(lambda x: pathogens.index(x), output_data["level_1"]))
    print(output_data)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=output_data, x="x", y="y", hue="correlation", size="-log10(p)", legend="brief", ax=ax)

    matplotlib.pyplot.xlabel("Pathogens")
    matplotlib.pyplot.ylabel("Pathogens")
    matplotlib.pyplot.xticks([])
    matplotlib.pyplot.yticks([])

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
