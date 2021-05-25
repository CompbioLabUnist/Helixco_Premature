"""
step01-5.py: Make manifest file for everything
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("metadata", help="Metadata XLSX file", type=str)
    parser.add_argument("--first", help="Input FASTQ.gz files", action="append", default=[])
    parser.add_argument("--second", help="Input FASTQ.gz files", action="append", default=[])
    parser.add_argument("--third", help="Input FASTQ.gz files", action="append", default=[])
    parser.add_argument("--stool", help="Input FASTQ.gz files", action="append", default=[])

    args = parser.parse_args()

    if not args.metadata.endswith(".xlsx"):
        raise ValueError("Metadata file must end with .xlsx!!")

    metadata = pandas.read_excel(args.metadata, dtype=str)

    args.first.sort()
    args.second.sort()
    args.third.sort()
    args.stool.sort()

    print("sample-id\tforward-absolute-filepath\treverse-absolute-filepath")

    f1 = list(filter(lambda x: x.endswith("1.fastq.gz"), args.first))
    f2 = list(filter(lambda x: x.endswith("2.fastq.gz"), args.first))
    assert len(f1) == len(f2)
    for fa, fb in zip(f1, f2):
        name = fa[:-11]
        name = name[name.rfind("/") + 1:].replace("_", "-").split("-")

        if len(name) == 3:
            name = "First-" + "-".join(name)
        elif (len(name) == 2) and (name[1] in ["B1", "B3", "B5"]):
            name = "First-" + name[0] + "-0-" + name[1]
        elif len(name) == 2:
            name = "First-M-" + "-".join(name)
        else:
            raise Exception("Something went wrong!!")

        print(name, fa, fb, sep="\t")

    f1 = list(filter(lambda x: x.endswith("1.fastq.gz"), args.second))
    f2 = list(filter(lambda x: x.endswith("2.fastq.gz"), args.second))
    assert len(f1) == len(f2)
    for fa, fb in zip(f1, f2):
        name = fa[:-11]
        name = name[name.rfind("/") + 1:].replace("_", "-").split("-")
        data = metadata.loc[(metadata["Sample"] == name[0]) & (metadata["ID"] == name[1]), :].to_numpy()[0]

        if data[1] == "Mother":
            name = "Second-" + data[0] + "-M-" + data[2]
        else:
            name = "Second-" + "-".join(data[:3])

        print(name, fa, fb, sep="\t")

    f1 = list(filter(lambda x: x.endswith("1.fastq.gz"), args.third))
    f2 = list(filter(lambda x: x.endswith("2.fastq.gz"), args.third))
    assert len(f1) == len(f2)
    for fa, fb in zip(f1, f2):
        name = fa[:-11]
        name = name[name.rfind("/") + 1:].replace("_", "-").split("-")
        data = metadata.loc[(metadata["Sample"] == name[0]) & (metadata["ID"] == name[1]), :].to_numpy()[0]

        if data[1] == "Mother":
            name = "Third-" + data[0] + "-M-" + data[2]
        else:
            name = "Third-" + "-".join(data[:3])

        print(name, fa, fb, sep="\t")

    f1 = list(filter(lambda x: x.endswith("R1_001.fastq.gz"), args.stool))
    f2 = list(filter(lambda x: x.endswith("R2_001.fastq.gz"), args.stool))
    assert len(f1) == len(f2)
    for fa, fb in zip(f1, f2):
        name = fa[:-19]
        name = name[name.rfind("/") + 1:].replace("_", "-").split("-")[3]
        data = metadata.loc[(metadata["Sample"].isin(["S1", "S3", "S5"])) & (metadata["ID"] == name), :].to_numpy()[0]
        name = "Stool-" + "-".join(data[:3])
        print(name, fa, fb, sep="\t")
