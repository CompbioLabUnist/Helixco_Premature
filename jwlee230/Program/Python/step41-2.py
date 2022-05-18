"""
step41-2.py: Heatmap Plot for correlation between Bacteria & Clinical data with site separation
"""
import argparse
import tarfile
import multiprocessing
import matplotlib
import matplotlib.pyplot
import pandas
import scipy.stats
import seaborn
import tqdm
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
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pearson", help="Pearson r", action="store_true", default=False)
    group.add_argument("--spearman", help="Spearman r", action="store_true", default=False)
    group.add_argument("--kendall", help="kendall tau", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tar.gz"):
        raise ValueError("Input file must end with .tar.gz!!")
    elif not args.metadata.endswith(".tsv"):
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

    merged_data = pandas.concat([input_data, metadata], axis="columns", verify_integrity=True, join="inner")
    print(merged_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()

    for site in tqdm.tqdm(step00.selected_long_sites):
        output_data = pandas.DataFrame(index=sorted(step00.numeric_columns), columns=taxa)
        input_data = merged_data.loc[(metadata["Site"] == site)]

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

        g = seaborn.clustermap(data=output_data, figsize=(48, 18), row_cluster=False, col_cluster=True, xticklabels=False, yticklabels=True, square=False, cmap="coolwarm", vmin=-1, vmax=1, cbar_pos=(-0.04, 0.2, 0.02, 0.6))

        g.ax_heatmap.set_xlabel("{0} Pathogens".format(len(taxa)))

        figures.append("{0}.pdf".format(site))
        g.savefig(figures[-1])

    with tarfile.open(args.output, "w") as tar:
        for file in tqdm.tqdm(figures):
            tar.add(file, arcname=file)
