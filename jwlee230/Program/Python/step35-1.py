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


def pearson(a: str, b: str) -> typing.Union[float, None]:
    correlation, p = scipy.stats.pearsonr(input_data[a], input_data[b])
    if p < 0.01:
        return correlation
    else:
        return None


def spearman(a: str, b: str) -> typing.Union[float, None]:
    correlation, p = scipy.stats.spearmanr(input_data[a], input_data[b])
    if p < 0.01:
        return correlation
    else:
        return None


def kendall(a: str, b: str) -> typing.Union[float, None]:
    correlation, p = scipy.stats.kendalltau(input_data[a], input_data[b])
    if p < 0.01:
        return correlation
    else:
        return None


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

    input_data = step00.read_pickle(args.input).T
    input_data.columns = list(map(step00.simplified_taxonomy, list(input_data.columns)))
    input_data = input_data.iloc[1:, :]
    print(input_data)

    output_data = pandas.DataFrame(data=numpy.zeros((len(input_data.columns), len(input_data.columns))), index=list(input_data.columns), columns=list(input_data.columns), dtype=float)

    filtered_columns = list(filter(lambda x: len(x.split(";")) > 5, list(output_data.columns)))
    output_data = output_data.loc[filtered_columns, filtered_columns]

    with multiprocessing.Pool(args.cpus) as pool:
        for a in list(output_data.index):
            if args.pearson:
                output_data.loc[a, :] = pool.starmap(pearson, [(a, b) for b in list(output_data.index)])
            elif args.spearman:
                output_data.loc[a, :] = pool.starmap(spearman, [(a, b) for b in list(output_data.index)])
            elif args.kendall:
                output_data.loc[a, :] = pool.starmap(kendall, [(a, b) for b in list(output_data.index)])
            else:
                raise Exception("Something went wrong!!")

    output_data = output_data.dropna(axis="index", thresh=len(output_data.index) // 10).fillna(0)
    output_data = output_data.loc[output_data.index, output_data.index]
    print(output_data)

    g = seaborn.clustermap(data=output_data, figsize=(32, 32), row_cluster=True, col_cluster=True, cbar_pos=None, xticklabels=False, yticklabels=False, square=False, cmap="bwr", vmin=-1, center=0, vmax=1)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")

    g.savefig(args.output)
