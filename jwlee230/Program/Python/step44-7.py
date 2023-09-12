"""
step44-7.py: Venn diagram for differentially abundant taxa
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import upsetplot
import tqdm
import step00

ratio_threshold = 2
p_threshold = 0.05

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--annotation", help="Annotated name(s)", type=str, nargs="+")

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    elif len(args.input) != len(args.annotation):
        raise ValueError("Length of annotation must match length of input!!")

    input_dict = dict()
    for file, annot in tqdm.tqdm(zip(args.input, args.annotation)):
        input_data = pandas.read_csv(file, sep="\t", index_col=0)
        input_dict[f"{annot}-Up"] = set(input_data.loc[(input_data["log2FoldChange"] > numpy.log2(ratio_threshold)) & (input_data["padj"] < p_threshold)].index)
        input_dict[f"{annot}-Down"] = set(input_data.loc[(input_data["log2FoldChange"] < (-1 * numpy.log2(ratio_threshold))) & (input_data["padj"] < p_threshold)].index)

    print("Union:", len(set.union(*input_dict.values())))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig = matplotlib.pyplot.figure(figsize=(10 * len(input_dict) + 10, 24))

    if set.union(*input_dict.values()):
        upsetplot.plot(upsetplot.from_contents(input_dict), fig=fig, show_counts="%d", show_percentages=True, element_size=None)

    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output, bbox_inches="tight")
    matplotlib.pyplot.close(fig)
