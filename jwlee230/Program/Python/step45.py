"""
step45.py: alpha-diversity
"""
import argparse
import io
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import skbio
import skbio.diversity
import statannotations.Annotator
import tqdm
import step00

metadata = pandas.DataFrame()


def draw_all(meta: str, alpha: str) -> str:
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    order = sorted(set(metadata[meta]))

    try:
        seaborn.violinplot(data=metadata, x=meta, y=alpha, ax=ax, order=order, linewidth=10, cut=1)
    except TypeError:
        matplotlib.pyplot.close(fig)
        return ""

    try:
        statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, r=2)), data=metadata, x=meta, y=alpha, order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
    except ValueError:
        pass

    matplotlib.pyplot.title("All")
    matplotlib.pyplot.ylabel(alpha.replace("_", " "))
    matplotlib.pyplot.tight_layout()

    fig_name = f"All+{meta.replace(' ', '_')}+{alpha}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


def draw_alpha(site: str, meta: str, alpha: str) -> str:
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    output_data = metadata.loc[(metadata["Site"] == site)]
    order = sorted(set(output_data[meta]))

    try:
        seaborn.violinplot(data=output_data, x=meta, y=alpha, ax=ax, order=order, linewidth=10, cut=1)
    except TypeError:
        matplotlib.pyplot.close(fig)
        return ""

    try:
        statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, r=2)), data=output_data, x=meta, y=alpha, order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
    except ValueError:
        pass

    matplotlib.pyplot.title(site)
    matplotlib.pyplot.ylabel(alpha.replace("_", " "))
    matplotlib.pyplot.tight_layout()

    fig_name = f"{site}+{meta.replace(' ', '_')}+{alpha}.pdf"
    fig.savefig(fig_name)
    matplotlib.pyplot.close(fig)

    return fig_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("tree", type=str, help="Tree NWk file")
    parser.add_argument("output", type=str, help="Output TAR f ile")
    parser.add_argument("--cpus", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if args.cpus < 1:
        raise ValueError("CPUS must be greater than zero")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID")
    del input_data["taxonomy"]
    input_data = input_data.T
    print(input_data)

    with open(args.tree, "r") as f:
        tree = skbio.TreeNode.read(io.StringIO(f.readline()))

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[sorted(set(input_data.index) & set(metadata.index)), :].replace(to_replace=-1, value=None)
    metadata_columns = sorted(set(metadata.columns) - step00.numeric_columns - {"Mother", "Too much weight gain", "Site"})
    print(metadata)
    print(metadata_columns)

    alphas = ["faith_pd", "observed_otus", "pielou_e", "shannon"]
    for alpha in tqdm.tqdm(alphas):
        if alpha.endswith("_ci"):
            continue

        try:
            if alpha == "faith_pd":
                metadata[alpha] = skbio.diversity.alpha_diversity(alpha, input_data, ids=list(input_data.index), otu_ids=list(input_data.columns), tree=tree)
            else:
                metadata[alpha] = skbio.diversity.alpha_diversity(alpha, input_data, ids=list(input_data.index))
        except (TypeError, KeyError, ValueError) as e:
            print(alpha, e)

    print(metadata)
    print(alphas)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    with multiprocessing.Pool(args.cpus) as pool:
        figures = pool.starmap(draw_all, itertools.product(metadata_columns, alphas))
        figures += pool.starmap(draw_alpha, itertools.product(step00.selected_long_sites, metadata_columns, alphas))
        figures = list(filter(None, figures))

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
