"""
step39.py: Beta-diversity heatmap plot
"""
import argparse
import itertools
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot
import pandas
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    sites = sorted(set(map(lambda x: x.split("-")[-1], list(input_data.index))))
    colorings = dict(zip(sites, itertools.cycle(matplotlib.colors.TABLEAU_COLORS)))
    site_colors = list(map(lambda x: colorings[x.split("-")[-1]], list(input_data.index)))
    print(colorings)

    g = seaborn.clustermap(data=input_data, figsize=(32, 32), row_cluster=True, col_cluster=True, row_colors=site_colors, col_colors=site_colors, xticklabels=False, yticklabels=False, cmap="Reds_r")

    matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=colorings[site]) for site in sites], sites, title="Sites", bbox_to_anchor=(1, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure, loc="best")

    g.savefig(args.output)
