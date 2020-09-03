"""
step10: select only useful features
"""
import argparse
import typing
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


def read_ancom_data(path: str) -> typing.List[str]:
    """
    read_ancom_data: Read & return data from ANCOM
    """
    if not path.endswith(".tsv"):
        raise ValueError("Input is not TSV file")

    data = pandas.read_csv(path, sep="\t")
    data = data.loc[(data["Reject null hypothesis"])]

    return list(data[list(data.columns)[0]])


def select_data(raw_data: pandas.DataFrame, significance: typing.List[str]) -> pandas.DataFrame:
    """
    select_data: select data from raw data by ANCOM
    """
    significance = list(map(step00.consistency_taxonomy, significance))

    return raw_data[list(filter(lambda x: step00.consistency_taxonomy(x[0]) in significance, list(raw_data.columns)))]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("raw", help="Raw TSV file", type=str, nargs=1)
    parser.add_argument("ancom", help="ANCOM TSV file", type=str, nargs=1)
    parser.add_argument("output", help="Output TAR.gz file", type=str, nargs=1)

    args = parser.parse_args()

    raw_data = read_raw_data(args.raw[0])
    ancom_data = read_ancom_data(args.ancom[0])
    data = select_data(raw_data, ancom_data)

    step00.make_pickle(args.output[0], data)
