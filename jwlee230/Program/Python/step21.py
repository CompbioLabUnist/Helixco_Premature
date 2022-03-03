"""
step21.py: Read & Clearify raw TSV into pandas
"""
import argparse
import pandas
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("output", help="Output TAR.gz file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")

    raw_data = pandas.read_csv(args.input, sep="\t")
    print(raw_data)

    output_data = raw_data.groupby(["taxonomy"]).sum()
    print(output_data)

    step00.make_pickle(args.output, output_data)
