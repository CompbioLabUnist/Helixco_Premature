"""
step51.py: QIIME2 - Batch corrected
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tsv file", type=str)
    parser.add_argument("output", help="Output tsv file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    if not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .tsv!!")

    data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=["#OTU ID", "taxonomy"])
    print(data)
    data.to_csv(args.output, sep="\t", index_label=("Feature ID", "Taxon"))
