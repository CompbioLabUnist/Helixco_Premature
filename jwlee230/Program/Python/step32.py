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
import step00

distance_data = pandas.DataFrame()
data = pandas.DataFrame()


def draw(disease: str) -> str:
    print(disease)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(36, 36))
    seaborn.scatterplot(data=data, x="tSNE1", y="tSNE2", ax=ax, hue="Detail Premature", style=disease, hue_order=sorted(set(data["Detail Premature"])), markers={x: y for x, y in zip(sorted(set(data[disease])), itertools.cycle(step00.markers))}, s=40 ** 2)

    try:
        p_value = skbio.stats.distance.permanova(skbio.stats.distance.DistanceMatrix(distance_data), list(data[disease]))["p-value"]
    except ValueError:
        p_value = 1.0
    matplotlib.pyplot.title("{0} (p={1:.3f})".format(disease, p_value))

    file_name = "{0}.pdf".format(disease.replace(" ", "_"))
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

    distance_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(distance_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[list(distance_data.index), :].replace(to_replace=-1, value=None)
    diseases = list(metadata.columns)
    print(metadata)

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="exact", n_jobs=args.cpus, perplexity=50, n_iter=10 ** 5, verbose=1).fit_transform(distance_data), columns=["tSNE1", "tSNE2"])

    for column in list(tsne_data.columns):
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])
    tsne_data["index"] = list(distance_data.index)
    tsne_data.set_index(keys="index", inplace=True, verify_integrity=True)
    print(tsne_data)

    data = pandas.concat(objs=[tsne_data, metadata], axis="columns", verify_integrity=True)
    print(data)

    with multiprocessing.Pool(args.cpus) as pool:
        files = pool.map(draw, diseases)

    with tarfile.open(args.output, "w") as tar:
        for f in files:
            print("Compressing:", f)
            tar.add(f, arcname=f)
