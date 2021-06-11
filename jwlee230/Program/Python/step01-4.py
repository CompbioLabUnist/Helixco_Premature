"""
step01-4.py: make manifest file from stool data
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("metadata", help="Metadata XLSX file", type=str)
    parser.add_argument("--stool", help="Input FASTQ.gz files", action="append", default=[])

    args = parser.parse_args()

    if not args.metadata.endswith(".xlsx"):
        raise ValueError("Metadata file must end with .xlsx!!")

    args.stool.sort()

    metadata = pandas.read_excel(args.metadata, dtype=str)

    f1 = sorted(list(filter(lambda x: x.endswith("R1_001.fastq.gz"), args.stool)))
    f2 = sorted(list(filter(lambda x: x.endswith("R2_001.fastq.gz"), args.stool)))
    assert len(f1) == len(f2)

    print("sample-id\tforward-absolute-filepath\treverse-absolute-filepath")
    assert len(f1) == len(f2)
    for fa, fb in zip(f1, f2):
        name = fa[:-19]
        name = name[name.rfind("/") + 1:].replace("_", "-").split("-")[3]
        data = metadata.loc[(metadata["Sample"].isin(["S1", "S3", "S5"])) & (metadata["ID"] == name), :].to_numpy()[0]
        name = "Stool-" + "-".join(data[:3])
        print(name, fa, fb, sep="\t")
