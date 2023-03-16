"""
step85.py: Antibiotics effects
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
column = "Mother Antibiotics"


def draw(metric: str, site: str) -> str:
    selected_distance_data = distance_data[metric].loc[(metadata["Site"] == site), (metadata["Site"] == site)].copy()
    selected_tsne_data = tsne_data[metric].loc[(tsne_data[metric]["Site"] == site)].copy(deep=True)

    for c in ["tSNE1", "tSNE2"]:
        selected_tsne_data[c] = sklearn.preprocessing.scale(selected_tsne_data[c])

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    seaborn.scatterplot(data=selected_tsne_data, x="tSNE1", y="tSNE2", ax=ax, hue="Premature", hue_order=["PTB", "Normal"], s=40 ** 2, palette=step00.PTB_two_colors, style=column, style_order=[True, False], legend="brief")

    try:
        p_value = skbio.stats.distance.permanova(skbio.stats.distance.DistanceMatrix(selected_distance_data), list(selected_tsne_data[column]), permutations=step00.small)["p-value"]
    except ValueError:
        p_value = 1.0

    matplotlib.pyplot.title(f"{column} (p={p_value:.3f})")
    matplotlib.pyplot.tight_layout()

    file_name = f"Antibiotics+{metric}+{site}.pdf"
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


def draw_no(metric: str, site: str) -> str:
    selected_distance_data = distance_data[metric].loc[(metadata["Site"] == site), (metadata["Site"] == site)].copy()
    selected_tsne_data = tsne_data[metric].loc[(tsne_data[metric]["Site"] == site)].copy(deep=True)

    for c in ["tSNE1", "tSNE2"]:
        selected_tsne_data[c] = sklearn.preprocessing.scale(selected_tsne_data[c])

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    seaborn.scatterplot(data=selected_tsne_data, x="tSNE1", y="tSNE2", ax=ax, hue="Premature", hue_order=["PTB", "Normal"], s=40 ** 2, palette=step00.PTB_two_colors, legend="brief")

    for item, color in step00.PTB_two_colors.items():
        step00.confidence_ellipse(selected_tsne_data.loc[(selected_tsne_data["Premature"] == item), "tSNE1"], selected_tsne_data.loc[(selected_tsne_data["Premature"] == item), "tSNE2"], ax, facecolor=color, alpha=0.3)

    try:
        p_value = skbio.stats.distance.permanova(skbio.stats.distance.DistanceMatrix(selected_distance_data), list(selected_tsne_data["Premature"]), permutations=step00.small)["p-value"]
    except ValueError:
        p_value = 1.0

    matplotlib.pyplot.title(f"PTB/Normal (p={p_value:.3f})")
    matplotlib.pyplot.tight_layout()

    file_name = f"All+{metric}+{site}.pdf"
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output PDF file")
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

    input_data = pandas.read_csv(args.input, sep="\t", index_col="#OTU ID", skiprows=1)
    input_data["taxonomy"] = list(map(lambda x: step00.consistency_taxonomy(x) if (step00.filtering_taxonomy(x)) else "Unclassified", input_data["taxonomy"]))
    input_data = input_data.groupby("taxonomy").sum().T
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col="#SampleID")
    metadata = metadata.loc[sorted(set(input_data.index) & set(metadata.index)), :]
    print(metadata)

    input_data = input_data.loc[sorted(set(input_data.index) & set(metadata.index)), :]
    print(input_data)

    for metric in tqdm.tqdm(step00.pdist_list):
        try:
            distance_data[metric] = skbio.diversity.beta_diversity(metric, input_data.to_numpy(), list(input_data.index)).to_data_frame()
        except Exception:
            continue

        tsne_data[metric] = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="exact", n_jobs=args.cpus, perplexity=50).fit_transform(distance_data[metric]), columns=["tSNE1", "tSNE2"])

        tsne_data[metric]["index"] = list(distance_data[metric].index)
        tsne_data[metric].set_index(keys="index", inplace=True, verify_integrity=True)

        for c in [column, "Site", "Premature"]:
            tsne_data[metric][c] = list(map(lambda x: metadata.loc[x, c], list(tsne_data[metric].index)))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    files = list()
    with multiprocessing.Pool(args.cpus) as pool:
        files += pool.starmap(draw, itertools.product(distance_data.keys(), step00.selected_long_sites))
        files += pool.starmap(draw_no, itertools.product(distance_data.keys(), step00.selected_long_sites))

    files = list(filter(None, files))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(files):
            tar.add(f, arcname=f)
