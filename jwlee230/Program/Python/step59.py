"""
step59.py: Violin Plots for Taxonomy Distribution
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00

metadata = pandas.DataFrame()
drawing_data = pandas.DataFrame()


def query(index, column):
    return metadata.loc[index, column]


def draw_all(meta: str) -> str:
    order = sorted(set(drawing_data[meta]))

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    try:
        seaborn.violinplot(data=drawing_data, x=meta, y="Taxonomy", order=order, linewidth=5, cut=1, ax=ax)
    except TypeError:
        matplotlib.pyplot.close(fig)
        return ""

    try:
        statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, r=2)), data=drawing_data, x=meta, y="Taxonomy", order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
    except ValueError:
        pass

    matplotlib.pyplot.title("All")
    matplotlib.pyplot.ylabel("Taxonomy Abundances")
    matplotlib.pyplot.tight_layout()

    fig_name = f"All+{meta.replace(' ', '_')}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


def draw(meta: str, site: str) -> str:
    tmp_data = drawing_data.loc[(drawing_data["Site"] == site)]
    order = sorted(set(tmp_data[meta]))

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    try:
        seaborn.violinplot(data=tmp_data, x=meta, y="Taxonomy", order=order, linewidth=5, cut=1, ax=ax)
    except TypeError:
        matplotlib.pyplot.close(fig)
        return ""

    try:
        statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, r=2)), data=tmp_data, x=meta, y="Taxonomy", order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
    except ValueError:
        pass

    matplotlib.pyplot.title(site)
    matplotlib.pyplot.ylabel("Taxonomy Abundances")
    matplotlib.pyplot.tight_layout()

    fig_name = f"{site}+{meta.replace(' ', '_')}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output TAR file")
    parser.add_argument("--cpus", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if args.cpus < 1:
        raise ValueError("CPUS must be greater than zero")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID")
    del input_data["taxonomy"]
    input_data = input_data.T
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[sorted(set(input_data.index) & set(metadata.index)), :].replace(to_replace=-1, value=None)
    metadata_columns = sorted(set(metadata.columns) - step00.numeric_columns - {"Mother", "Too much weight gain", "Site"})
    print(metadata)
    print(metadata_columns)

    drawing_data = pandas.DataFrame(data=[(index, input_data.loc[index, column]) for index, column in itertools.product(list(input_data.index), list(input_data.columns))], columns=["Sample", "Taxonomy"])
    print(drawing_data)

    with multiprocessing.Pool(args.cpus) as pool:
        for column in tqdm.tqdm(metadata.columns):
            drawing_data[column] = pool.starmap(query, [(index, column) for index in drawing_data["Sample"]])
    print(drawing_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = pool.map(draw_all, metadata_columns)
        figures += pool.starmap(draw, itertools.product(metadata_columns, step00.selected_long_sites))
    figures = list(filter(None, figures))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
