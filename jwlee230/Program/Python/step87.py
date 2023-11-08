"""
step87.py: Select DAT for FTB-DAT & PROM-DAT
"""
import argparse
import numpy
import pandas

ratio_threshold = 2
p_threshold = 0.05


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("FTB", help="DAT TSV file", type=str,)
    parser.add_argument("PROM", help="DAT TSV file", type=str,)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.FTB.endswith(".tsv"):
        raise ValueError("FTB file must end with.TSV!!")
    elif not args.PROM.endswith(".tsv"):
        raise ValueError("PROM file must end with.TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end with.TSV!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    FTB_data = pandas.read_csv(args.FTB, sep="\t", index_col=0)
    FTB_data = FTB_data.loc[((FTB_data["log2FoldChange"] > numpy.log2(ratio_threshold)) | (FTB_data["log2FoldChange"] < numpy.log2(1 / ratio_threshold))) & (FTB_data["padj"] < p_threshold)]
    print(FTB_data)

    PROM_data = pandas.read_csv(args.PROM, sep="\t", index_col=0)
    PROM_data = PROM_data.loc[((PROM_data["log2FoldChange"] > numpy.log2(ratio_threshold)) | (PROM_data["log2FoldChange"] < numpy.log2(1 / ratio_threshold))) & (PROM_data["padj"] < p_threshold)]
    print(PROM_data)

    selected_DAT = sorted(set(FTB_data.index) - set(PROM_data.index))
    print("DAT:", len(selected_DAT))

    input_data = input_data.loc[selected_DAT, :]
    print(input_data)
    input_data.to_csv(args.output, sep="\t")
