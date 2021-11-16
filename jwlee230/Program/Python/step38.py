"""
step38.py: Clustermap plot with sites
"""
import argparse
import itertools
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tar.gz"):
        raise ValueError("Input file must end with .tar.gz!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUS must be a positive integer!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = step00.read_pickle(args.input)
    input_data.index = list(map(step00.simplified_taxonomy, list(input_data.index)))
    input_data.sort_index(inplace=True)
    input_data = input_data.groupby(input_data.index).sum().T
    print(input_data)

    for index in tqdm.tqdm(list(input_data.index)):
        input_data.loc[index, :] = input_data.loc[index, :] / sum(input_data.loc[index, :])

    sites = sorted(set(map(lambda x: x.split("-")[-1], list(input_data.index))))
    colorings = dict(zip(sites, itertools.cycle(matplotlib.colors.TABLEAU_COLORS)))
    site_colors = list(map(lambda x: colorings[x.split("-")[-1]], list(input_data.index)))
    print(colorings)

    g = seaborn.clustermap(data=input_data, figsize=(32, 18), row_cluster=True, col_cluster=True, row_colors=site_colors, xticklabels=False, yticklabels=False, cmap="Reds", vmin=0, vmax=1)
    g.ax_heatmap.set_xlabel("Taxonomy")
    g.ax_heatmap.set_ylabel("Sample")

    g.savefig(args.output)
