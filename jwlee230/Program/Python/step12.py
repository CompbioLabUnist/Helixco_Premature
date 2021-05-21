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

    parser.add_argument("input", type=str, help="Input TAR.gz file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output PDF file")

    args = parser.parse_args()

    if not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .tsv")

    data: pandas.DataFrame = step00.read_pickle(args.input)
    info_data = pandas.read_csv(args.metadata, sep="\t", skiprows=[1])

    hue_style = args.output.split("/")[-1].split(".")[0]
    if "_" in hue_style:
        hue, style = hue_style.split("_")
    else:
        hue = style = hue_style

    print(hue, style)

    data[hue] = info_data[hue]
    data[style] = info_data[style]

    print(data)

    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=data, x="tSNE1", y="tSNE2", ax=ax, hue=hue, style=style, legend="brief", s=1000, hue_order=sorted(data[hue]))

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
