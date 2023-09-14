"""
step85.py: DAT violin plot
"""
import argparse
import tarfile
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

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("DAT", type=str, help="DAT TSV file", nargs="+")
    parser.add_argument("output", type=str, help="Output TAR file")
    parser.add_argument("--site", type=str, help="Site to plot", required=True)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif list(filter(lambda x: not x.endswith(".tsv"), args.DAT)):
        raise ValueError("DAT file must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID")
    input_data = input_data.groupby("taxonomy").sum().T
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str, index_col="#SampleID").replace(to_replace=-1, value=None).dropna(axis="columns", how="all")
    metadata = metadata.loc[(metadata["Site"] == args.site)]
    print(metadata)

    input_data = pandas.concat([input_data, metadata], axis="columns", join="inner", verify_integrity=True)
    print(input_data)

    DAT_list = list()
    for DAT_file in tqdm.tqdm(args.DAT):
        DAT_data = pandas.read_csv(DAT_file, sep="\t", index_col=0)
        DAT_list.append(set(DAT_data.loc[((DAT_data["log2FoldChange"] < numpy.log2(1 / ratio_threshold)) | (DAT_data["log2FoldChange"] > numpy.log2(ratio_threshold))) & (DAT_data["padj"] < p_threshold)].index))
    DAT = set.intersection(*DAT_list)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for taxon in tqdm.tqdm(DAT):
        simple_name = step00.simplified_taxonomy(taxon)

        fig, axs = matplotlib.pyplot.subplots(figsize=(36, 18), ncols=2, sharey=True)

        PTB_data = input_data.loc[(input_data["Premature"] == "PTB")]
        seaborn.violinplot(data=PTB_data, x="PROM", y=taxon, order=["False", "True"], palette=["tab:blue", "tab:red"], linewidth=5, ax=axs[0])
        axs[0].set_xlabel("PROM in PTB")
        axs[0].set_ylabel(simple_name)

        seaborn.violinplot(data=input_data, x="Premature", y=taxon, order=["Normal", "PTB"], palette=[step00.PTB_two_colors["Normal"], step00.PTB_two_colors["PTB"]], linewidth=5, ax=axs[1])
        axs[1].set_xlabel("PTB")
        axs[1].set_ylabel(simple_name)

        matplotlib.pyplot.tight_layout()

        figures.append(f"{simple_name}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for figure in tqdm.tqdm(figures):
            tar.add(figure, arcname=figure)
