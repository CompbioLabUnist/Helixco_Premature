"""
step55.py: Clustermap for Taxonomy Abundance
"""
import argparse
import pandas
import matplotlib
import matplotlib.pyplot
import seaborn
import step00
if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("output", type=str, help="Output PDF file")

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=["#OTU ID"]).groupby("taxonomy").sum().T
    input_data = input_data.loc[:, list(filter(step00.filtering_taxonomy, list(input_data.columns)))]
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    g = seaborn.clustermap(data=input_data, figsize=(32, 18), row_cluster=True, col_cluster=True, xticklabels=False, yticklabels=False, square=False, cmap="Reds", robust=True)

    g.ax_heatmap.set_xlabel(f"{input_data.shape[1]} bacteria")
    g.ax_heatmap.set_ylabel(f"{input_data.shape[0]} samples")

    g.savefig(args.output)
    matplotlib.pyplot.close(g.fig)
