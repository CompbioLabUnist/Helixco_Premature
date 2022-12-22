"""
step63-1.py: select metadata for metagenomeSeq
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--site", help="Select site", type=str, required=True)

    comparing = parser.add_mutually_exclusive_group()
    comparing.add_argument("--EL", help="Early & Late", action="store_true", default=False)
    comparing.add_argument("--EF", help="Early & Normal", action="store_true", default=False)
    comparing.add_argument("--LF", help="Late & Normal", action="store_true", default=False)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end with .TSV!!")

    metadata = pandas.read_csv(args.input, sep="\t", skiprows=[1], index_col=0).sort_values("Detail Premature")
    print(metadata)

    metadata = metadata.loc[(metadata["Site"] == args.site)]
    print(metadata)

    if args.EL:
        metadata = metadata.loc[(metadata["Detail Premature"].isin({"Early PTB", "Late PTB"}))]
        metadata["Comparing"] = metadata["Detail Premature"]
    elif args.EF:
        metadata = metadata.loc[(metadata["Detail Premature"].isin({"Early PTB", "Normal"}))]
        metadata["Comparing"] = metadata["Detail Premature"]
    elif args.LF:
        metadata = metadata.loc[(metadata["Detail Premature"].isin({"Late PTB", "Normal"}))]
        metadata["Comparing"] = metadata["Detail Premature"]
    else:
        metadata["Comparing"] = metadata["Simple Premature"]
    print(metadata)

    metadata.to_csv(args.output, sep="\t")
