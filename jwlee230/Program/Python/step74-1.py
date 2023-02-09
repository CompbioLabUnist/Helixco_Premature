"""
step74-1.py: draw relplot
"""
import argparse
import matplotlib
import matplotlib.pyplot
import seaborn
import pandas
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Train TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output PDF file")

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    print(input_data)

    taxa = list(map(step00.simplified_taxonomy, list(input_data)))
    input_data.columns = taxa
    input_data = input_data.loc[:, sorted(taxa)]

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    input_data = input_data.loc[sorted(set(input_data.index) & set(metadata.index))]
    print(input_data)

    corr_mat = input_data.corr().stack().reset_index(name="correlation")
    print(corr_mat)

    g = seaborn.relplot(data=corr_mat, x="level_0", y="level_1", hue="correlation", size="correlation", palette="coolwarm", hue_norm=(-1, 1), height=48, aspect=1, sizes=(250, 1250), size_norm=(-1, 1))
    g.set(xlabel="", ylabel="", aspect="equal")
    g.despine(left=True, bottom=True)
    g.ax.margins(0.02)
    for label in g.ax.get_xticklabels():
        label.set_rotation(90)
    g.tight_layout()
    g.fig.savefig(args.output)
