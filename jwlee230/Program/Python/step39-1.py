"""
step39-1.py: Beta-diversity heatmap plot for Macrogen data
"""
import argparse
import itertools
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot
import seaborn
import skbio.diversity
import skbio.tree
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--beta", choices=skbio.diversity.get_beta_diversity_metrics(), required=True)

    args = parser.parse_args()

    if not args.input.endswith(".tar.gz"):
        raise ValueError("Input file must end with .tar.gz!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = step00.read_pickle(args.input)
    input_data.index = list(map(step00.simplified_taxonomy, list(input_data.index)))
    input_data = input_data.iloc[1:, 1:].T
    print(input_data)

    sites = sorted(set(map(lambda x: x.split("-")[-1], list(input_data.index))))
    colorings = dict(zip(sites, itertools.cycle(matplotlib.colors.TABLEAU_COLORS)))
    site_colors = list(map(lambda x: colorings[x.split("-")[-1]], list(input_data.index)))
    print(colorings)

    tree = skbio.tree.TreeNode.from_taxonomy([(x, step00.simplified_taxonomy(x).split(";")) for x in list(input_data.columns)])

    for e in tqdm.tqdm(tree.traverse()):
        if e.is_root():
            continue
        e.length = len(e.name.split(";"))

    distance_data = skbio.diversity.beta_diversity(args.beta, input_data.to_numpy(), list(input_data.index), otu_ids=list(input_data.columns), tree=tree).to_data_frame()
    print(distance_data)

    g = seaborn.clustermap(data=distance_data, figsize=(32, 32), row_cluster=True, col_cluster=True, row_colors=site_colors, col_colors=site_colors, xticklabels=False, yticklabels=False, cmap="Reds_r")

    matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=colorings[site]) for site in sites], sites, title="Sites", bbox_to_anchor=(1, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure, loc="best")

    g.savefig(args.output)
