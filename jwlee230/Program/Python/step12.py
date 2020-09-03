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
    parser.add_argument("output", type=str, nargs=1, help="Output PNG file")

    args = parser.parse_args()

    if not args.output[0].endswith(".png"):
        raise ValueError("Output file must end with .png")

    data: pandas.DataFrame = step00.read_pickle(args.input[0])

    seaborn.set(context="poster", style="whitegrid")
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=data, x="TSNE1", y="TSNE2", ax=ax, legend=False)

    fig.savefig(args.output[0])
    matplotlib.pyplot.close(fig)
