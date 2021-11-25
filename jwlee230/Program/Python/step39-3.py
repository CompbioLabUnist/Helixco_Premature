"""
step39-3.py: Beta-diversity heatmap plot wth site for Macrogen data
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot
import pandas
import seaborn
import skbio.diversity
import skbio.tree
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)
    parser.add_argument("--beta", choices=skbio.diversity.get_beta_diversity_metrics(), required=True)

    args = parser.parse_args()

    if not args.input.endswith(".tar.gz"):
        raise ValueError("Input file must end with .tar.gz!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")

    input_data = step00.read_pickle(args.input)
    input_data.index = list(map(step00.simplified_taxonomy, list(input_data.index)))
    input_data = input_data.iloc[1:, 1:].T
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col=0)
    print(metadata)
    print(sorted(metadata.columns))

    tree = skbio.tree.TreeNode.from_taxonomy([(x, step00.simplified_taxonomy(x).split(";")) for x in list(input_data.columns)])

    for e in tqdm.tqdm(tree.traverse()):
        if e.is_root():
            continue
        e.length = len(e.name.split(";"))

    distance_data = skbio.diversity.beta_diversity(args.beta, input_data.to_numpy(), list(input_data.index), otu_ids=list(input_data.columns), tree=tree).to_data_frame()
    print(distance_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig_files = list()
    for site in tqdm.tqdm(set(metadata["Site"])):
        drawing_data = distance_data.loc[(metadata["Site"] == site), (metadata["Site"] == site)]

        colors = pandas.DataFrame(index=input_data.index)
        colors["Detail Premature"] = list(map(lambda x: step00.PTB_colors[x], metadata["Detail Premature"]))
        for disease in sorted(["Gestational Diabetes", "Too much weight gain", "Obesity", "Hypertension", "PROM", "Neonate Antibiotics", "Mother Antibiotics", "Mother Steroid"]):
            colors[disease] = list(map(lambda x: "k" if metadata.loc[x, disease] else "w", list(input_data.index)))

        g = seaborn.clustermap(data=drawing_data, figsize=(64, 32), row_cluster=True, col_cluster=True, row_colors=colors, col_colors=colors["Detail Premature"], xticklabels=False, yticklabels=False, cmap="Reds_r")

        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=step00.PTB_colors[PTB]) for PTB in step00.detailed_PTB], step00.detailed_PTB, title="Detailed PTB", bbox_to_anchor=(1, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure, loc="best")

        fig_files.append("{0}.pdf".format(site))
        g.savefig(fig_files[-1])

    with tarfile.open(args.output, "w") as tar:
        for file in tqdm.tqdm(fig_files):
            tar.add(file, arcname=file)
