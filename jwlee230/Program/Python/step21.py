"""
step21.py: Read & Clearify raw TSV into pandas
"""
import argparse
import pandas


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    raw_data = pandas.read_csv(args.input, sep="\t", skiprows=1)
    print(raw_data)

    output_data = raw_data.groupby(["taxonomy"]).mean()
    print(output_data)

    output_data.to_csv(args.output, sep="\t")
