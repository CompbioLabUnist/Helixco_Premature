"""
step39-2.py: Beta-diversity heatmap plot with site
"""
import argparse
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot
import pandas
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata must end with .tsv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col=0)
    print(metadata)
    print(sorted(metadata.columns))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig_files = list()
    for site in tqdm.tqdm(set(metadata["Site"])):
        drawing_data = input_data.loc[(metadata["Site"] == site), (metadata["Site"] == site)]

        colors = pandas.DataFrame(index=drawing_data.index)
        colors["Detail Premature"] = list(map(lambda x: step00.PTB_colors[metadata.loc[x, "Detail Premature"]], list(drawing_data.index)))
        for disease in sorted(["Gestational Diabetes", "Too much weight gain", "Obesity", "Hypertension", "PROM", "Neonate Antibiotics", "Mother Antibiotics", "Mother Steroid"]):
            colors[disease] = list(map(lambda x: "k" if metadata.loc[x, disease] else "w", list(drawing_data.index)))

        g = seaborn.clustermap(data=drawing_data, figsize=(64, 32), row_cluster=True, col_cluster=True, row_colors=colors, col_colors=colors["Detail Premature"], xticklabels=False, yticklabels=False, cmap="Reds_r")

        matplotlib.pyplot.legend([matplotlib.patches.Patch(facecolor=step00.PTB_colors[PTB]) for PTB in step00.detailed_PTB], step00.detailed_PTB, title="Detailed PTB", bbox_to_anchor=(1, 1), bbox_transform=matplotlib.pyplot.gcf().transFigure, loc="best")

        fig_files.append("{0}.pdf".format(site))
        g.savefig(fig_files[-1])

    with tarfile.open(args.output, "w") as tar:
        for file in tqdm.tqdm(fig_files):
            tar.add(file, arcname=file)
