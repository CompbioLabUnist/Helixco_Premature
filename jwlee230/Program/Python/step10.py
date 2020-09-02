"""
step10: select only useful features
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

    parser.add_argument("raw", help="Raw TSV file", type=str, nargs=1)
    parser.add_argument("ancom", help="ANCOM TSV file", type=str, nargs=1)

    args = parser.parse_args()

    data = read_raw_data(args.raw)
    print(data)
