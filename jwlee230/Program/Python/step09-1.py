"""
step09-1.py: make log
"""
import argparse
import numpy
import pandas


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("output", type=str, help="Output TSV file")

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID")
    input_data.iloc[:, :-1] = numpy.log(input_data.iloc[:, :-1] + 1)
    print(input_data)

    with open(args.output, "w") as f:
        f.write("# Constructed from TSV file\n")
        f.write(input_data.to_csv(sep="\t"))
