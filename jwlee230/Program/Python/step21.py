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

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=0)
    input_data["taxonomy"] = list(map(lambda x: x.replace(".", "-"), input_data["taxonomy"]))
    print(input_data)

    output_data = input_data.groupby(["taxonomy"]).sum()
    output_data.columns = list(map(lambda x: x.replace(".", "-"), output_data.columns))
    print(output_data)

    step00.make_pickle(args.output, output_data)
