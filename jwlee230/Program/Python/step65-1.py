"""
step65-1.py: Venn diagram for differentially abundant taxa
"""
import argparse
import numpy
import pandas
import tqdm

ratio_threshold = 2
p_threshold = 0.05


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--annotation", help="Annotated name(s)", type=str, nargs="+")

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end with .TSV!!")
    elif len(args.input) != len(args.annotation):
        raise ValueError("Length of annotation must match length of input!!")

    input_dict = dict()
    DAT_set = set()
    for file, annot in tqdm.tqdm(zip(args.input, args.annotation)):
        input_data = pandas.read_csv(file, sep="\t", index_col=0)

        up_data = input_data.loc[(input_data["log2FoldChange"] > numpy.log2(ratio_threshold)) & (input_data["padj"] < p_threshold)]
        input_dict[annot + "-Up"] = set(up_data.index)
        DAT_set |= set(up_data.index)

        down_data = input_data.loc[(input_data["log2FoldChange"] < numpy.log2(1 / ratio_threshold)) & (input_data["padj"] < p_threshold)]
        input_dict[annot + "-Down"] = set(down_data.index)
        DAT_set |= set(down_data.index)

    output_data = pandas.DataFrame(index=sorted(DAT_set), columns=args.annotation, dtype=str)
    output_data.loc[:, :] = ""
    for taxon in tqdm.tqdm(DAT_set):
        for annot in args.annotation:
            if taxon in input_dict[annot + "-Up"]:
                output_data.loc[taxon, annot] = "Up"
            elif taxon in input_dict[annot + "-Down"]:
                output_data.loc[taxon, annot] = "Down"

    print(output_data)
    output_data.to_csv(args.output, sep="\t")
