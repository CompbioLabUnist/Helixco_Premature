"""
step21: Read & Clearify raw TSV into pandas
"""
import argparse
import pandas
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str, nargs=1)
    parser.add_argument("output", help="Output TAR.gz file", type=str, nargs=1)

    args = parser.parse_args()

    raw_data = pandas.read_csv(args.input[0], sep="\t", skiprows=1)
    raw_data.set_index(inplace=True, keys=["taxonomy", "#OTU ID"], verify_integrity=True)
    raw_data = raw_data.T

    data = pandas.DataFrame()
    taxonomy_list = sorted(list(set(map(lambda x: x[0], raw_data.columns))))
    for taxonomy in taxonomy_list:
        data[taxonomy] = raw_data[list(filter(lambda x: x[0] == taxonomy, raw_data.columns))].sum(axis=1)

    step00.make_pickle(args.output[0], data)
