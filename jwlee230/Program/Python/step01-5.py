"""
step01-5.py: Make manifest file for everything
"""
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--first", help="Input FASTQ.gz files", action="append", default=[])
    parser.add_argument("--second", help="Input FASTQ.gz files", action="append", default=[])
    parser.add_argument("--third", help="Input FASTQ.gz files", action="append", default=[])
    parser.add_argument("--stool", help="Input FASTQ.gz files", action="append", default=[])

    args = parser.parse_args()

    args.first.sort()
    args.second.sort()
    args.third.sort()
    args.stool.sort()

    print("sample-id\tforward-absolute-filepath\treverse-absolute-filepath")

    f1 = sorted(list(filter(lambda x: x.endswith("1.fastq.gz"), args.first)))
    f2 = sorted(list(filter(lambda x: x.endswith("2.fastq.gz"), args.first)))
    assert len(f1) == len(f2)
    for fa, fb in zip(f1, f2):
        name = fa[:-11]
        if fb.startswith(name):
            print(name[name.rfind("/") + 1:].replace("_", "-"), fa, fb, sep="\t")
        else:
            raise ValueError(fa, fb)

    f1 = sorted(list(filter(lambda x: x.endswith("1.fastq.gz"), args.second)))
    f2 = sorted(list(filter(lambda x: x.endswith("2.fastq.gz"), args.second)))
    assert len(f1) == len(f2)
    for fa, fb in zip(f1, f2):
        name = fa[:-11]
        if fb.startswith(name):
            print(name[name.rfind("/") + 1:].replace("_", "-"), fa, fb, sep="\t")
        else:
            raise ValueError(fa, fb)

    f1 = sorted(list(filter(lambda x: x.endswith("R1_001.fastq.gz"), args.stool)))
    f2 = sorted(list(filter(lambda x: x.endswith("R2_001.fastq.gz"), args.stool)))
    assert len(f1) == len(f2)
    for fa, fb in zip(f1, f2):
        name = fa[:-19]
        if fb.startswith(name):
            print(name[name.rfind("/") + 1:].replace("_", "-"), fa, fb, sep="\t")
        else:
            raise ValueError(fa, fb)
