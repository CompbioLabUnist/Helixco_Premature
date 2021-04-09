"""
step10-2.py: select everything and make pandas
"""
import argparse
import pandas
import step00


def read_raw_data(path: str) -> pandas.DataFrame:
    """
    read_raw_data: Read & clearify data from raw TSV
    """
    if not path.endswith(".tsv"):
        raise ValueError("Input is not TSV file")

    data = pandas.read_csv(path, sep="\t", skiprows=1)
    data.set_index(inplace=True, keys=["taxonomy", "#OTU ID"], verify_integrity=True)
    data = data.T

    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Raw TSV file", type=str)
    parser.add_argument("output", help="Output TAR.gz file", type=str)

    args = parser.parse_args()

    step00.make_pickle(args.output, read_raw_data(args.input))
