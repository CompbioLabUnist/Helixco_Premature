"""
step12.py: draw basic t-SNE
"""
import argparse
import pandas
import matplotlib
import matplotlib.pyplot
import seaborn
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, nargs=1, help="Input TAR.gz file")
    parser.add_argument("metadata", type=str, nargs=1, help="Metadata TSV file")
    parser.add_argument("output", type=str, nargs=1, help="Output PNG file")

    args = parser.parse_args()

    if not args.output[0].endswith(".png"):
        raise ValueError("Output file must end with .png")
    elif not args.metadata[0].endswith(".tsv"):
        raise ValueError("Metadata file must end with .tsv")

    data: pandas.DataFrame = step00.read_pickle(args.input[0])
    info_data = pandas.read_csv(args.metadata[0], sep="\t", skiprows=[1])

    hue_style = args.output[0].split("/")[-1].split(".")[0]
    if "_" in hue_style:
        hue, style = hue_style.split("_")
    else:
        hue = style = hue_style

    data[hue] = info_data[hue]
    data[style] = info_data[style]

    seaborn.set(context="poster", style="whitegrid")
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=data, x="TSNE1", y="TSNE2", ax=ax, hue=hue, style=style, legend="brief")

    fig.savefig(args.output[0])
    matplotlib.pyplot.close(fig)
