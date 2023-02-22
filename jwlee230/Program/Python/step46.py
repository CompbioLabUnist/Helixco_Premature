"""
step46.py: beta-diversity
"""
import argparse
import io
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import skbio
import skbio.diversity
import skbio.stats
import sklearn.manifold
import sklearn.preprocessing
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("tree", type=str, help="Tree NWk file")
    parser.add_argument("output", type=str, help="Output TAR file")
    parser.add_argument("--cpus", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if args.cpus < 1:
        raise ValueError("CPUS must be greater than zero")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.tree.endswith(".nwk"):
        raise ValueError("Tree must end with .NWK!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID")
    del input_data["taxonomy"]
    input_data = input_data.T
    print(input_data)

    with open(args.tree, "r") as f:
        tree = skbio.TreeNode.read(io.StringIO(f.readline()))

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[sorted(set(input_data.index) & set(metadata.index)), :].replace(to_replace=-1, value=None).dropna(axis="columns")
    metadata_columns = sorted(set(metadata.columns) - step00.numeric_columns - {"Mother", "Too much weight gain", "Site"})
    print(metadata)
    print(metadata_columns)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for beta in skbio.diversity.get_beta_diversity_metrics():
        beta_diversity = skbio.diversity.beta_diversity(beta, input_data, ids=list(input_data.index), otu_ids=list(input_data.columns), tree=tree)

        tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(init="pca", random_state=42, verbose=1, method="exact", n_jobs=args.cpus).fit_transform(beta_diversity.to_data_frame()), columns=["tSNE1", "tSNE2"])
        for column in tqdm.tqdm(list(tsne_data.columns)):
            metadata[column] = sklearn.preprocessing.scale(tsne_data[column])
        print(metadata)

        for meta in tqdm.tqdm(metadata_columns):
            p = skbio.stats.distance.permanova(beta_diversity, grouping=metadata, column=meta)["p-value"]
            order = sorted(set(metadata[meta]))

            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

            seaborn.scatterplot(data=metadata, x="tSNE1", y="tSNE2", hue=meta, hue_order=order, ax=ax, s=1000)

            matplotlib.pyplot.title(f"All (PERMANOVA p={p:.3f})")
            matplotlib.pyplot.tight_layout()

            fig_name = f"All+{meta.replace(' ', '_')}+{beta}.pdf"
            fig.savefig(fig_name)
            figures.append(fig_name)
            matplotlib.pyplot.close(fig)

        beta_diversity_df = beta_diversity.to_data_frame()
        for site, meta in tqdm.tqdm(list(itertools.product(step00.selected_long_sites, metadata_columns))):
            tmp_meta = metadata.loc[(metadata["Site"] == site)]
            order = sorted(set(tmp_meta[meta]))
            beta_diversity = skbio.stats.distance.DistanceMatrix(beta_diversity_df.loc[list(tmp_meta.index), list(tmp_meta.index)], ids=list(tmp_meta.index))

            try:
                p = skbio.stats.distance.permanova(beta_diversity, grouping=tmp_meta, column=meta)["p-value"]
            except ValueError:
                p = 1.0

            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

            seaborn.scatterplot(data=tmp_meta, x="tSNE1", y="tSNE2", hue=meta, hue_order=order, ax=ax, s=1000)

            matplotlib.pyplot.title(f"{site} (PERMANOVA p={p:.3f}")
            matplotlib.pyplot.tight_layout()

            fig_name = f"{site}+{meta.replace(' ', '_')}+{beta}.pdf"
            fig.savefig(fig_name)
            figures.append(fig_name)
            matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
