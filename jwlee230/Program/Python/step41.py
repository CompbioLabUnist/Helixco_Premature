"""
step41.py: Heatmap Plot for correlation between Bacteria & Clinical data
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import step00

input_data = pandas.DataFrame()


def pearson(clinical: str, taxonomy: str) -> float:
    drawing_data = input_data[[taxonomy, clinical]]
    drawing_data = drawing_data.loc[(drawing_data[clinical] != -1)]
    return scipy.stats.pearsonr(drawing_data[taxonomy], drawing_data[clinical])[0]


def spearman(clinical: str, taxonomy: str) -> float:
    drawing_data = input_data[[taxonomy, clinical]]
    drawing_data = drawing_data.loc[(drawing_data[clinical] != -1)]
    return scipy.stats.spearmanr(drawing_data[taxonomy], drawing_data[clinical])[0]


def kendall(clinical: str, taxonomy: str) -> float:
    drawing_data = input_data[[taxonomy, clinical]]
    drawing_data = drawing_data.loc[(drawing_data[clinical] != -1)]
    return scipy.stats.kendalltau(drawing_data[taxonomy], drawing_data[clinical])[0]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pearson", help="Pearson r", action="store_true", default=False)
    group.add_argument("--spearman", help="Spearman r", action="store_true", default=False)
    group.add_argument("--kendall", help="kendall tau", action="store_true", default=False)

    args = parser.parse_args()

    if not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata must end with .tsv!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    elif args.cpus < 0:
        raise ValueError("CPUs must be positive!!")

    input_data = step00.read_pickle(args.input)
    input_data.index = list(map(step00.simplified_taxonomy, list(input_data.index)))
    input_data.sort_index(inplace=True)
    input_data = input_data.groupby(input_data.index).sum().T
    taxa = list(filter(lambda x: x.count(";") > 5, list(input_data.columns)))
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col=0)
    print(metadata)
    print(sorted(metadata.columns))

    input_data = pandas.concat([input_data, metadata], axis="columns", verify_integrity=True, join="inner")
    print(input_data)

    output_data = pandas.DataFrame(index=sorted(step00.numeric_columns), columns=taxa)
    print(output_data)

    with multiprocessing.Pool(args.cpus) as pool:
        for index in list(output_data.index):
            if args.pearson:
                output_data.loc[index, :] = pool.starmap(pearson, [(index, column) for column in list(output_data.columns)])
            elif args.spearman:
                output_data.loc[index, :] = pool.starmap(spearman, [(index, column) for column in list(output_data.columns)])
            elif args.kendall:
                output_data.loc[index, :] = pool.starmap(kendall, [(index, column) for column in list(output_data.columns)])
            else:
                raise Exception("Something went wrong!!")

    output_data.fillna(0, inplace=True)
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(48, 18))

    seaborn.heatmap(data=output_data, vmin=-1, vmax=1, cmap="coolwarm", xticklabels=False, yticklabels=True, ax=ax)

    matplotlib.pyplot.xlabel("Pathogens")

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
