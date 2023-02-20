"""
step44-4.py: draw box plots
"""
import argparse
import itertools
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

ratio_threshold = 2
p_threshold = 0.05

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("DAT", help="Output DAT file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col=0)
    print(metadata)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    input_data["Premature"] = list(map(lambda x: metadata.loc[x, "Premature"], input_data.index))
    print(input_data)

    DAT_data = pandas.read_csv(args.DAT, sep="\t", index_col=0)
    print(DAT_data)

    DAT_data = DAT_data.loc[((DAT_data["log2FoldChange"] > numpy.log2(ratio_threshold)) | (DAT_data["log2FoldChange"] < numpy.log2(1 / ratio_threshold))) & (DAT_data["padj"] < p_threshold)].sort_values("padj")
    print(DAT_data)

    taxa_list = list(DAT_data.index)
    raw_violin_data: typing.List[typing.Tuple[str, str, float]] = list()
    for taxon in tqdm.tqdm(taxa_list):
        raw_violin_data.extend(zip(itertools.cycle([taxon]), input_data["Premature"], input_data[taxon]))
    violin_data = pandas.DataFrame(raw_violin_data, columns=["Feature", "Premature", "Value"])
    print(violin_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(len(taxa_list) * 1.5, 24))
    seaborn.boxplot(data=violin_data, x="Feature", y="Value", hue="Premature", order=taxa_list, hue_order=("PTB", "Normal"), palette=step00.PTB_two_colors, showfliers=False, ax=ax)
    matplotlib.pyplot.xlabel("")
    matplotlib.pyplot.xticks(range(len(taxa_list)), list(map(step00.simplified_taxonomy, taxa_list)), rotation=90)
    matplotlib.pyplot.tight_layout()
    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
