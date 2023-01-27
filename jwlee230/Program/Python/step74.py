"""
step74.py: draw pairplot
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

    taxa = list(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    input_data = pandas.concat([input_data, metadata], axis="columns", join="inner", verify_integrity=True)
    print(input_data)

    g = seaborn.pairplot(data=input_data, hue="Detail Premature", hue_order=step00.detailed_PTB, palette=step00.PTB_colors, x_vars=taxa, y_vars=taxa, height=4, aspect=1, kind="reg", diag_kind="kde")
    g.tight_layout()
    g.fig.savefig(args.output)
