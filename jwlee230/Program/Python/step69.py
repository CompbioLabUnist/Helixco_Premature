"""
step69.py: Select DAT of DESeq2 for classification
"""
import argparse
import numpy
import pandas
import tqdm
import step00

ratio_threshold = 2
p_threshold = 0.05


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("DAT", help="DAT TSV file(s)", type=str, nargs="+")
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif list(filter(lambda x: not x.endswith(".tsv"), args.DAT)):
        raise ValueError("DAT file must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end with .TSV!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    DAT = list(input_data.index)
    for file in tqdm.tqdm(args.DAT):
        DAT_data = pandas.read_csv(file, sep="\t", index_col=0)
        DAT_data = DAT_data.loc[list(map(step00.filtering_taxonomy, list(DAT_data.index))), :]
        DAT_data = DAT_data.loc[((DAT_data["log2FoldChange"] > numpy.log2(ratio_threshold)) | (DAT_data["log2FoldChange"] < numpy.log2(1 / ratio_threshold))) & (DAT_data["padj"] < p_threshold)]
        DAT = list(filter(lambda x: x in DAT, list(DAT_data.index)))

    print("Union:", len(DAT))

    input_data = input_data.loc[DAT, :]
    print(input_data)
    input_data.to_csv(args.output, sep="\t")
