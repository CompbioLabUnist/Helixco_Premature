"""
step32.py: Draw beta-diversity scatter plots
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import skbio.stats.distance
import sklearn.manifold
import sklearn.preprocessing
import tqdm
import step00

distance_data = pandas.DataFrame()
data = pandas.DataFrame()


def draw(disease: str, site: str) -> str:
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    tmp_data = data.loc[(data["Site"] == site)]
    seaborn.scatterplot(data=tmp_data, x="tSNE1", y="tSNE2", ax=ax, hue=disease, hue_order=sorted(set(tmp_data[disease])), s=40 ** 2, edgecolor=None)

    try:
        tmp_distance_data = distance_data.loc[(data["Site"] == site), (data["Site"] == site)]
        p_value = 1.0
        if not tmp_distance_data.empty:
            p_value = skbio.stats.distance.permanova(skbio.stats.distance.DistanceMatrix(tmp_distance_data), list(tmp_data[disease]), permutations=step00.small)["p-value"]
    except ValueError:
        p_value = 1.0
    matplotlib.pyplot.title("{0} (p={1:.3f})".format(disease, p_value))

    file_name = "{0}+{1}.pdf".format(site, disease.replace(" ", "_"))
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


def draw_all(disease: str) -> str:
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    seaborn.scatterplot(data=data, x="tSNE1", y="tSNE2", ax=ax, hue=disease, hue_order=sorted(set(data[disease])), s=40 ** 2, edgecolor=None)

    try:
        p_value = skbio.stats.distance.permanova(skbio.stats.distance.DistanceMatrix(distance_data), list(data[disease]), permutations=step00.small)["p-value"]
    except ValueError:
        p_value = 1.0
    matplotlib.pyplot.title("{0} (p={1:.3f})".format(disease, p_value))

    file_name = "All+{0}.pdf".format(disease.replace(" ", "_"))
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output TAR file")
    parser.add_argument("--cpus", type=int, default=1, help="CPU to use")

    data_group = parser.add_mutually_exclusive_group()
    data_group.add_argument("--first", help="Select First data", action="store_true", default=False)
    data_group.add_argument("--second", help="Select Second+Third data", action="store_true", default=False)

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

    distance_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    distance_data = distance_data.iloc[:min(distance_data.shape), :min(distance_data.shape)]

    if args.first:
        distance_data = distance_data.loc[list(filter(lambda x: x.startswith("First"), list(distance_data.index))), list(filter(lambda x: x.startswith("First"), list(distance_data.columns)))]
    elif args.second:
        distance_data = distance_data.loc[list(filter(lambda x: x.startswith("Second") or x.startswith("Third"), list(distance_data.index))), list(filter(lambda x: x.startswith("Second") or x.startswith("Third"), list(distance_data.columns)))]

    print(distance_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[sorted(set(distance_data.index) & set(metadata.index)), sorted(set(metadata.columns) - step00.numeric_columns)].replace(to_replace=-1, value=None)
    diseases = set(metadata.columns) - step00.numeric_columns - {"Mother", "Neonate", "Site"}
    print(metadata)
    print(sorted(diseases))

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="exact", n_jobs=args.cpus, perplexity=50, n_iter=step00.big, verbose=1).fit_transform(distance_data), columns=["tSNE1", "tSNE2"])

    for column in list(tsne_data.columns):
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])
    tsne_data["index"] = list(distance_data.index)
    tsne_data.set_index(keys="index", inplace=True, verify_integrity=True)
    print(tsne_data)

    data = pandas.concat(objs=[tsne_data, metadata], axis="columns", join="inner", verify_integrity=True)
    print(data)

    with multiprocessing.Pool(args.cpus) as pool:
        files = pool.starmap(draw, itertools.product(diseases, step00.selected_long_sites))
        files += pool.map(draw_all, diseases)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(files):
            tar.add(f, arcname=f)
