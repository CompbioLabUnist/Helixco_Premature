"""
step42.py: Violin plot for taxonomy with clinical data
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import statannotations.Annotator
import tqdm
import step00

data = pandas.DataFrame()
threshold = 0.01


def draw(taxo: str, disease: str, site: str) -> str:
    drawing_data = data.loc[(data["Site"] == site)]
    order = sorted(set(drawing_data[disease]))

    if len(order) < 1:
        return ""

    for a in order:
        min_value = min(drawing_data.loc[(drawing_data[disease] == a), taxo])
        median_value = numpy.median(drawing_data.loc[(drawing_data[disease] == a), taxo])
        max_value = max(drawing_data.loc[(drawing_data[disease] == a), taxo])

        if (min_value == median_value) or (median_value == max_value):
            return ""

    for a, b in list(itertools.combinations(order, 2)):
        stat, p = scipy.stats.mannwhitneyu(drawing_data.loc[(drawing_data[disease] == a), taxo], drawing_data.loc[(drawing_data[disease] == b), taxo])
        if p < threshold:
            break
    else:
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.violinplot(data=drawing_data, x=disease, y=taxo, order=order, inner="box", linewidth=5, cut=1, ax=ax)
    statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, 2)), data=drawing_data, x=disease, y=taxo, order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
    matplotlib.pyplot.scatter(x=range(len(order)), y=[numpy.mean(drawing_data.loc[(drawing_data[disease] == d), taxo]) for d in order], marker="*", c="white", s=400, zorder=10)

    matplotlib.pyplot.title(site)
    matplotlib.pyplot.ylabel(step00.consistency_taxonomy(taxo, number=1))
    matplotlib.pyplot.tight_layout()

    file_name = "{2}+{1}+{0}.pdf".format(step00.consistency_taxonomy(taxo).replace(";", "_"), disease.replace(" ", "_"), site)
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


def draw_all(taxo: str, disease: str) -> str:
    order = sorted(set(data[disease]))

    if len(order) < 1:
        return ""

    for a in order:
        min_value = min(data.loc[(data[disease] == a), taxo])
        median_value = numpy.median(data.loc[(data[disease] == a), taxo])
        max_value = max(data.loc[(data[disease] == a), taxo])

        if (min_value == median_value) or (median_value == max_value):
            return ""

    for a, b in list(itertools.combinations(order, 2)):
        stat, p = scipy.stats.mannwhitneyu(data.loc[(data[disease] == a), taxo], data.loc[(data[disease] == b), taxo])
        if p < threshold:
            break
    else:
        return ""

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.violinplot(data=data, x=disease, y=taxo, order=order, inner="box", linewidth=5, cut=1, ax=ax)
    statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, 2)), data=data, x=disease, y=taxo, order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
    matplotlib.pyplot.scatter(x=range(len(order)), y=[numpy.mean(data.loc[(data[disease] == d), taxo]) for d in order], marker="*", c="white", s=400, zorder=10)

    matplotlib.pyplot.ylabel(step00.consistency_taxonomy(taxo, number=1))
    matplotlib.pyplot.tight_layout()

    file_name = "{2}+{1}+{0}.pdf".format(step00.consistency_taxonomy(taxo).replace(";", "_"), disease.replace(" ", "_"), "All")
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TAR.gz file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output TAR file")
    parser.add_argument("--cpus", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if not args.input.endswith(".tar.gz"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUs must be greater than zero!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=0).groupby("taxonomy").sum().T
    taxonomy_list = list(input_data.columns)
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[list(set(input_data.index) & set(metadata.index)), sorted(set(metadata.columns) - step00.numeric_columns)].replace(to_replace=-1, value=None)
    diseases = set(metadata.columns) - step00.numeric_columns - {"Mother", "Neonate", "Site"}
    print(metadata)
    print(sorted(diseases))

    data = pandas.concat(objs=[input_data, metadata], axis="columns", join="inner", verify_integrity=True)
    print(data)

    with multiprocessing.Pool(processes=args.cpus) as pool:
        figures = list(filter(None, pool.starmap(draw, itertools.product(taxonomy_list, diseases, step00.selected_long_sites))))
        figures += list(filter(None, pool.starmap(draw_all, itertools.product(taxonomy_list, diseases))))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
