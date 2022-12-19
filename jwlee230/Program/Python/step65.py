"""
step65.py: Venn diagram for differentially abundant taxa
"""
import argparse
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import tqdm
import venn
import step00

ratio_threshold = 2
p_threshold = 0.05


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--annotation", help="Annotated name(s)", type=str, nargs="+")

    updown = parser.add_mutually_exclusive_group(required=True)
    updown.add_argument("--up", action="store_true", default=False)
    updown.add_argument("--down", action="store_true", default=False)

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
        if args.up:
            input_data = input_data.loc[(input_data["log2FoldChange"] > numpy.log2(ratio_threshold)) & (input_data["padj"] < p_threshold)]
        elif args.down:
            input_data = input_data.loc[(input_data["log2FoldChange"] < numpy.log2(1 / ratio_threshold)) & (input_data["padj"] < p_threshold)]
        else:
            raise Exception("Something went wrong!!")
        input_dict[annot] = set(input_data.index)

    print("Union:", len(set.union(*input_dict.values())))

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(10, 10))

    venn.venn(input_dict, fmt=step00.venn_format, ax=ax)

    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
