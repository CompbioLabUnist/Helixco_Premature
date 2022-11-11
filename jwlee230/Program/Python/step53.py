"""
step53.py: Abundance Distribution
"""
import argparse
import itertools
import multiprocessing
import pandas
import matplotlib
import matplotlib.pyplot
import seaborn
import step00

data = pandas.DataFrame()
metadata = pandas.DataFrame()


def run(index, column):
    return (metadata.loc[index, "Site"], data.loc[index, column])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output PDF file")
    parser.add_argument("--cpus", type=int, default=1, help="Number of CPUs to use")

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("Number of CPUs must be positive!!")

    data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=["#OTU ID"]).T.iloc[:-1, :]
    print(data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    with multiprocessing.Pool(args.cpus) as pool:
        output_data = pandas.DataFrame(pool.starmap(run, itertools.product(list(data.index), list(data.columns))), columns=["Site", "Abundance"])
    output_data = output_data.loc[(output_data["Abundance"] > 0)]
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

    seaborn.histplot(data=output_data, x="Abundance", hue="Site", stat="probability", kde=True, multiple="stack", ax=ax)

    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
