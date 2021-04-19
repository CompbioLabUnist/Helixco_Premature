"""
step22.py: draw t-SNE for brief
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
    parser.add_argument("output", type=str, help="Output PDF file")

    args = parser.parse_args()

    if not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")

    data: pandas.DataFrame = step00.read_pickle(args.input)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=data, x="tSNE1", y="tSNE2", ax=ax, legend=None)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
