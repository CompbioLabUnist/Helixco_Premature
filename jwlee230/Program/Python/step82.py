"""
step82.py: Correlation within DAT
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
    parser.add_argument("output", type=str, help="Output PDF file")

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    print(input_data)

    taxa = list(map(step00.simplified_taxonomy, list(input_data.columns)))
    input_data.columns = taxa
    print(input_data)

    corr_data = input_data.corr().stack().reset_index(name="correlation")
    corr_data.loc[(abs(corr_data["correlation"]) < 0.5), "correlation"] = 0
    print(corr_data)

    g = seaborn.relplot(data=corr_data, x="level_0", y="level_1", hue="correlation", size="correlation", palette="vlag", hue_norm=(-1, 1), edgecolor=".7", height=48, sizes=(800, 4000), size_norm=(-1, 1))

    g.set(xlabel="", ylabel="", aspect="equal")
    g.tick_params(labelsize="xx-small")
    g.despine(left=True, bottom=True)

    for label in g.ax.get_xticklabels():
        label.set_rotation(90)
    for artist in g.legend.legendHandles:
        artist.set_edgecolor(".7")

    g.tight_layout()
    g.fig.savefig(args.output)
