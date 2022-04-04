"""
ste31.py: draw violin plot with ANCOM
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("ancom", help="ANCOM ancom.tsv file", type=str)
    parser.add_argument("data", help="Data tsv file", type=str)
    parser.add_argument("metadata", help="Metadata tsv file", type=str)
    parser.add_argument("output", help="Output TAR file", type=str)

    args = parser.parse_args()

    if not args.ancom.endswith(".tsv"):
        raise ValueError("ANCOM file must end with .TSV!!")

    ancom_data = pandas.read_csv(args.ancom, sep="\t", names=["id", "W", "Reject null hypothesis"], usecols=["id", "Reject null hypothesis"], header=0, index_col="id")
    print(ancom_data)

    input_data = pandas.read_csv(args.data, sep="\t", skiprows=1, index_col=0).groupby("taxonomy").sum().T
    print(input_data)

    print(sorted(ancom_data.index)[:5])
    print(sorted(input_data.columns)[:5])

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(how="all", axis="columns")
    print(metadata)
