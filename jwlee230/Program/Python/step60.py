"""
step60.py: Abundance/Proportion Distribution
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.colors
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output TAR file")

    args = parser.parse_args()

    if not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output must end with .TAR!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID").groupby("taxonomy").sum().T
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str, index_col="#SampleID").dropna(axis="columns", how="all")
    metadata = metadata.loc[sorted(set(input_data.index) & set(metadata.index)), :].replace(to_replace=-1, value=None).sort_values("Detail Gestational Week")
    print(metadata)

    input_data = input_data.loc[metadata.index, :]
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for s, site in step00.selected_sites_dict.items():
        tmp_data = input_data.loc[list(filter(lambda x: x.endswith(s), list(input_data.index))), :]
        taxa = sorted(list(tmp_data.columns), key=lambda x: sum(tmp_data[x]), reverse=True)
        tmp_data = tmp_data.loc[:, taxa]

        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

        for i, (color, taxon) in tqdm.tqdm(list(enumerate(zip(itertools.cycle(matplotlib.colors.XKCD_COLORS), taxa)))):
            if i < 5:
                try:
                    label = step00.simplified_taxonomy(taxon)
                except IndexError:
                    label = step00.consistency_taxonomy(taxon, 1)
                matplotlib.pyplot.bar(range(tmp_data.shape[0]), tmp_data.iloc[:, i], bottom=numpy.sum(tmp_data.iloc[:, :i], axis=1), color=color, linewidth=0, label=label)
            else:
                matplotlib.pyplot.bar(range(tmp_data.shape[0]), tmp_data.iloc[:, i], bottom=numpy.sum(tmp_data.iloc[:, :i], axis=1), color=color, linewidth=0)

        matplotlib.pyplot.xticks(range(tmp_data.shape[0]), list(map(lambda x: metadata.loc[x, "Detail Gestational Week"], list(tmp_data.index))), rotation="vertical", fontsize="xx-small")
        matplotlib.pyplot.xlabel(f"{tmp_data.shape[0]} {site} samples (order by GW)")
        matplotlib.pyplot.ylabel(f"{len(taxa)} bacteria")
        matplotlib.pyplot.legend()
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.tight_layout()

        figures.append(f"{site}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
