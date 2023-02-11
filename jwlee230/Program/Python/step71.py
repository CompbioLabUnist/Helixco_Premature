"""
step71.py: Volcano plot with differetially abundant taxa from metagenomeSeq
"""
import argparse
import itertools
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("metadata", help="Metadata TSV file", type=str)
    parser.add_argument("EL", help="EL data TSV file", type=str)
    parser.add_argument("EF", help="EF data TSV file", type=str)
    parser.add_argument("LF", help="LF data TSV file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("METADATA must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures: typing.List[str] = list()

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    taxa = list(input_data.columns)
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    input_data = input_data.loc[sorted(set(input_data.index) & set(metadata.index))]
    print(input_data)

    input_data["Detail Premature"] = list(map(lambda x: metadata.loc[x, "Detail Premature"], list(input_data.index)))
    print(input_data)

    EL_data = pandas.read_csv(args.EL, sep="\t", index_col=0)
    EF_data = pandas.read_csv(args.EF, sep="\t", index_col=0)
    LF_data = pandas.read_csv(args.LF, sep="\t", index_col=0)

    for taxon in tqdm.tqdm(taxa):
        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=input_data, x="Detail Premature", y=taxon, order=step00.detailed_PTB, linewidth=5, cut=1, ax=ax)

        try:
            statannotations.Annotator.Annotator(ax, list(itertools.combinations(step00.detailed_PTB, r=2)), data=input_data, x="Detail Premature", y=taxon, order=step00.detailed_PTB).configure(test=None, text_format="simple", loc="inside", comparisons_correction=None, verbose=0).set_pvalues_and_annotate([EL_data.loc[taxon, "pvalues"], EF_data.loc[taxon, "pvalues"], LF_data.loc[taxon, "pvalues"]])
        except ValueError:
            pass

        matplotlib.pyplot.scatter(x=range(len(step00.detailed_PTB)), y=[numpy.mean(input_data.loc[(input_data["Detail Premature"] == d), taxon]) for d in step00.detailed_PTB], marker="*", c="white", s=400, zorder=10)

        matplotlib.pyplot.xlabel("")
        matplotlib.pyplot.ylabel(step00.simplified_taxonomy(taxon))
        matplotlib.pyplot.tight_layout()

        figures.append(f"{step00.simplified_taxonomy(taxon)}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(figures):
            tar.add(f, arcname=f)
