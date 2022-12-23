"""
step63-2.py: select metadata for metagenomeSeq
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end with .TSV!!")

    data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=0).groupby("taxonomy").sum()
    print(data)
    data.to_csv(args.output, sep="\t")
