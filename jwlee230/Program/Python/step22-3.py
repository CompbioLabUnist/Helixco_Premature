"""
step22-3.py: draw t-SNE plot
"""
import argparse
import itertools
import tarfile
import pandas
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import seaborn
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("output", type=str, help="Output TAR file")

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .tar!!")

    tar_files = list()

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    columns = set(input_data.columns)
    columns -= {"tSNE1", "tSNE2"}
    print(sorted(columns))

    numeric_columns = step00.numeric_columns & columns
    categorical_columns = columns - step00.numeric_columns
    print(sorted(numeric_columns))
    print(sorted(categorical_columns))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    PTB_palette = dict(zip(step00.detailed_PTB, matplotlib.colors.TABLEAU_COLORS))

    for site in set(input_data["Site"]):
        data = input_data.loc[(input_data["Site"] == site)]

        for c in tqdm.tqdm(numeric_columns):
            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

            seaborn.scatterplot(data=data, x="tSNE1", y="tSNE2", ax=ax, hue="Detail Premature", size=c, legend="brief", hue_order=step00.detailed_PTB, palette=PTB_palette)

            for PTB, color in PTB_palette.items():
                step00.confidence_ellipse(data.loc[(data["Detail Premature"] == PTB), "tSNE1"], data.loc[(data["Detail Premature"] == PTB), "tSNE2"], ax, facecolor=color, alpha=0.3)

            matplotlib.pyplot.tight_layout()
            tar_files.append("{0}+{1}.pdf".format(site, c.replace(" ", "_")))
            fig.savefig(tar_files[-1])
            matplotlib.pyplot.close(fig)

        for c in tqdm.tqdm(categorical_columns):
            palette = dict(zip(sorted(set(data[c])), matplotlib.colors.XKCD_COLORS))

            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

            seaborn.scatterplot(data=data, x="tSNE1", y="tSNE2", ax=ax, hue=c, style=c, legend="brief", markers={x: y for x, y in zip(sorted(set(data[c])), itertools.cycle(step00.markers))}, s=40 ** 2, style_order=sorted(set(data[c])), palette=palette)

            for item, color in palette.items():
                step00.confidence_ellipse(data.loc[(data[c] == item), "tSNE1"], data.loc[(data[c] == item), "tSNE2"], ax, facecolor=color, alpha=0.3)

            matplotlib.pyplot.tight_layout()
            tar_files.append("{0}+{1}.pdf".format(site, c.replace(" ", "_")))
            fig.savefig(tar_files[-1])
            matplotlib.pyplot.close(fig)

    for c in tqdm.tqdm(numeric_columns):
        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.scatterplot(data=input_data, x="tSNE1", y="tSNE2", ax=ax, hue="Detail Premature", size=c, legend="brief", hue_order=step00.detailed_PTB)

        for PTB, color in PTB_palette.items():
            step00.confidence_ellipse(input_data.loc[(input_data["Detail Premature"] == PTB), "tSNE1"], input_data.loc[(input_data["Detail Premature"] == PTB), "tSNE2"], ax, facecolor=color, alpha=0.3)

        matplotlib.pyplot.tight_layout()
        tar_files.append("{0}+{1}.pdf".format("All", c.replace(" ", "_")))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    for c in tqdm.tqdm(categorical_columns):
        palette = dict(zip(sorted(set(input_data[c])), matplotlib.colors.XKCD_COLORS))

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.scatterplot(data=input_data, x="tSNE1", y="tSNE2", ax=ax, hue=c, style=c, legend="brief", markers={x: y for x, y in zip(sorted(set(input_data[c])), itertools.cycle(step00.markers))}, s=40 ** 2, style_order=sorted(set(input_data[c])))

        for item, color in palette.items():
            step00.confidence_ellipse(input_data.loc[(input_data[c] == item), "tSNE1"], input_data.loc[(input_data[c] == item), "tSNE2"], ax, facecolor=color, alpha=0.3)

        matplotlib.pyplot.tight_layout()
        tar_files.append("{0}+{1}.pdf".format("All", c.replace(" ", "_")))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for file_name in tar_files:
            tar.add(file_name, arcname=file_name)
