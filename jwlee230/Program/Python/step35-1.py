"""
step35.py: Draw cluster map plot within bacteria
"""
import argparse
import multiprocessing
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import step00

input_data = pandas.DataFrame()


def pearson(a: str, b: str, threshold: float) -> typing.Union[float, None]:
    correlation, p = scipy.stats.pearsonr(input_data[a], input_data[b])
    if p < threshold:
        return correlation
    else:
        return None


def spearman(a: str, b: str, threshold: float) -> typing.Union[float, None]:
    correlation, p = scipy.stats.spearmanr(input_data[a], input_data[b])
    if p < threshold:
        return correlation
    else:
        return None


def kendall(a: str, b: str, threshold: float) -> typing.Union[float, None]:
    correlation, p = scipy.stats.kendalltau(input_data[a], input_data[b])
    if p < threshold:
        return correlation
    else:
        return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("output", help="Output PNG file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)
    parser.add_argument("--p", help="P-value threshold", type=float, default=0.01)

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
    elif not (0 < args.p < 1):
        raise ValueError("P-value must be (0, 1)")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = step00.read_pickle(args.input).T
    input_data.columns = list(map(step00.simplified_taxonomy, list(input_data.columns)))
    input_data = input_data.iloc[1:, :]
    print(input_data)

    output_data = pandas.DataFrame(data=numpy.zeros((len(input_data.columns), len(input_data.columns))), index=list(input_data.columns), columns=list(input_data.columns), dtype=float)

    with multiprocessing.Pool(args.cpus) as pool:
        for a in list(output_data.index):
            if args.pearson:
                output_data.loc[a, :] = pool.starmap(pearson, [(a, b, args.p) for b in list(output_data.index)])
            elif args.spearman:
                output_data.loc[a, :] = pool.starmap(spearman, [(a, b, args.p) for b in list(output_data.index)])
            elif args.kendall:
                output_data.loc[a, :] = pool.starmap(kendall, [(a, b, args.p) for b in list(output_data.index)])
            else:
                raise Exception("Something went wrong!!")
    output_data.fillna(0, inplace=True)
    print(output_data)

    g = seaborn.clustermap(data=output_data, figsize=(64, 64), row_cluster=True, col_cluster=True, cbar=True, xticklabels=False, yticklabels=False, square=True, cmap="coolwarm", vmin=-1, center=0, vmax=1)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")

    g.savefig(args.output)
