"""
step21.py: Read & Clearify raw TSV into pandas
"""
import argparse
import pandas
import tqdm


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end with .TSV!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID")
    columns = list(input_data.columns)[:-1]
    for column in tqdm.tqdm(columns):
        input_data.loc[:, column] = input_data.loc[:, column] / sum(input_data.loc[:, column])
    print(input_data)

    with open(args.output, "w") as f:
        f.write("# Constructed from TSV file (Proportion)\n")

    input_data.to_csv(args.output, sep="\t", mode="a")
