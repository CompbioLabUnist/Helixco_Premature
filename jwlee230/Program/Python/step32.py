"""
step32.py: Draw beta-diversity scatter plots
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import pandas
import seaborn
import skbio.stats.distance
import skbio.diversity
import sklearn.manifold
import sklearn.preprocessing
import tqdm
import step00

distance_data: typing.Dict[str, pandas.DataFrame] = dict()
tsne_data: typing.Dict[str, pandas.DataFrame] = dict()
metadata = pandas.DataFrame()


def draw(metric: str, column: str, site: str) -> str:
    selected_distance_data = distance_data[metric].loc[(metadata["Site"] == site), (metadata["Site"] == site)]
    selected_tsne_data = tsne_data[metric].loc[(tsne_data[metric]["Site"] == site)]

    palette = dict(zip(sorted(set(selected_tsne_data[column])), matplotlib.colors.XKCD_COLORS))

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    seaborn.scatterplot(data=selected_tsne_data, x="tSNE1", y="tSNE2", ax=ax, hue=column, hue_order=sorted(set(selected_tsne_data[column])), s=40 ** 2, edgecolor=None, palette=palette)

    for item, color in palette.items():
        step00.confidence_ellipse(selected_tsne_data.loc[(selected_tsne_data[column] == item), "tSNE1"], selected_tsne_data.loc[(selected_tsne_data[column] == item), "tSNE2"], ax, facecolor=color, alpha=0.3)

    try:
        p_value = skbio.stats.distance.permanova(skbio.stats.distance.DistanceMatrix(selected_distance_data), list(selected_tsne_data[column]), permutations=step00.small)["p-value"]
    except ValueError:
        p_value = 1.0

    matplotlib.pyplot.title("{0} ({1}; p={2:.3f})".format(column, metric, p_value))
    matplotlib.pyplot.tight_layout()

    file_name = "{2}+{0}+{1}.pdf".format(column.replace(" ", "_"), metric, site)
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


def draw_all(metric: str, column: str) -> str:
    palette = dict(zip(sorted(set(tsne_data[metric][column])), matplotlib.colors.XKCD_COLORS))

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=tsne_data[metric], x="tSNE1", y="tSNE2", ax=ax, hue=column, hue_order=sorted(set(tsne_data[metric][column])), s=40 ** 2, edgecolor=None, palette=palette)

    for item, color in palette.items():
        step00.confidence_ellipse(tsne_data[metric].loc[(tsne_data[metric][column] == item), "tSNE1"], tsne_data[metric].loc[(tsne_data[metric][column] == item), "tSNE2"], ax, facecolor=color, alpha=0.3)

    try:
        p_value = skbio.stats.distance.permanova(skbio.stats.distance.DistanceMatrix(distance_data[metric]), list(tsne_data[metric][column]), permutations=step00.small)["p-value"]
    except ValueError:
        p_value = 1.0

    matplotlib.pyplot.title("{0} ({1}; p={2:.3f})".format(column, metric, p_value))
    matplotlib.pyplot.tight_layout()

    file_name = "All+{0}+{1}.pdf".format(column.replace(" ", "_"), metric)
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output TAR file")
    parser.add_argument("--cpus", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUS must be greater than zero!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col="#OTU ID", skiprows=1)
    del input_data["taxonomy"]
    input_data = input_data.T
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[sorted(set(input_data.index) & set(metadata.index)), sorted(set(metadata.columns) - step00.numeric_columns)].replace(to_replace=-1, value=None)
    columns = set(metadata.columns) - step00.numeric_columns - {"Mother", "Neonate", "Site"}
    columns = {"Data", "Site"}
    print(metadata)
    print(sorted(columns))

    for metric in tqdm.tqdm(step00.pdist_list):
        distance_data[metric] = skbio.diversity.beta_diversity(metric, input_data.to_numpy(), list(input_data.index)).to_data_frame()

        tsne_data[metric] = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="exact", n_jobs=args.cpus, perplexity=50).fit_transform(distance_data[metric]), columns=["tSNE1", "tSNE2"])

        for column in list(tsne_data[metric].columns):
            tsne_data[metric][column] = sklearn.preprocessing.scale(tsne_data[metric][column])

        tsne_data[metric]["index"] = list(distance_data[metric].index)
        tsne_data[metric].set_index(keys="index", inplace=True, verify_integrity=True)

        for column in columns:
            tsne_data[metric][column] = list(map(lambda x: metadata.loc[x, column], list(tsne_data[metric].index)))

    files = list()
    with multiprocessing.Pool(args.cpus) as pool:
        files += pool.starmap(draw, itertools.product(step00.pdist_list, columns, step00.selected_long_sites))
        files += pool.starmap(draw_all, itertools.product(step00.pdist_list, columns))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(files):
            tar.add(f, arcname=f)
