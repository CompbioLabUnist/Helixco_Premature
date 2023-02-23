"""
step75.py: Hierarchical clustering with abundance
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
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

    taxa = list(map(step00.simplified_taxonomy, list(input_data.columns)))
    input_data.columns = taxa

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    input_data = pandas.concat([input_data, metadata], axis="columns", join="inner", verify_integrity=True).sort_values("Detail Gestational Week")
    col_colors = list(map(lambda x: step00.PTB_two_colors[input_data.loc[x, "Premature"]], list(input_data.index)))
    input_data = numpy.log10(input_data[taxa].T + 1)
    print(input_data)

    g = seaborn.clustermap(data=input_data, figsize=(18, 32), xticklabels=False, yticklabels=True, row_cluster=True, col_cluster=False, col_colors=col_colors, cmap="YlOrRd", dendrogram_ratio=(0.2, 0.01), cbar_pos=(0.90, 0.85, 0.05, 0.13))
    g.fig.savefig(args.output)
