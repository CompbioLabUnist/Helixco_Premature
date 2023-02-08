"""
step61.py: Volcano plot with differetially abundant taxa from DESeq2
"""
import argparse
import itertools
import tarfile
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn
import statannotations.Annotator
import tqdm
import step00

ratio_threshold = 2
p_threshold = 0.05


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("DAT", help="DAT TSV file", type=str)
    parser.add_argument("coldata", help="coldata TSV file (not necessarily TSV)", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.DAT.endswith(".tsv"):
        raise ValueError("DAT must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    input_data = input_data.loc[:, list(filter(step00.select_taxonomy, list(input_data.columns)))]
    print(input_data)

    DAT_data = pandas.read_csv(args.DAT, sep="\t", index_col=0)
    DAT_data = DAT_data.loc[list(filter(step00.select_taxonomy, list(DAT_data.index))), :]
    DAT_data = DAT_data.loc[((DAT_data["log2FoldChange"] > numpy.log2(ratio_threshold)) | (DAT_data["log2FoldChange"] < -1 * numpy.log2(ratio_threshold))) & (DAT_data["padj"] < p_threshold)]
    print(DAT_data)

    coldata_data = pandas.read_csv(args.coldata, sep="\t", index_col=0)
    input_data["condition"] = coldata_data["condition"]
    order = sorted(set(coldata_data["condition"]))
    print(coldata_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    figures = list()
    for index in tqdm.tqdm(list(DAT_data.index)):
        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

        seaborn.violinplot(data=input_data, x="condition", y=index, order=order, linewidth=5, cut=1, ax=ax)

        try:
            statannotations.Annotator.Annotator(ax, list(itertools.combinations(order, r=2)), data=input_data, x="condition", y=index, order=order).configure(test="Mann-Whitney", text_format="simple", loc="inside", comparisons_correction=None, verbose=0).apply_and_annotate()
        except ValueError:
            pass

        matplotlib.pyplot.scatter(x=range(len(order)), y=[numpy.mean(drawing_data.loc[(drawing_data["condition"] == d), index]) for d in order], marker="*", c="white", s=400, zorder=10)

        matplotlib.pyplot.xlabel()
        matplotlib.pyplot.ylabel(step00.simplified_taxonomy(index))
        matplotlib.pyplot.tight_layout()

        figures.append(f"{step00.simplified_taxonomy(index)}.pdf")
        fig.savefig(figures[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(set(figures)):
            tar.add(f, arcname=f)
