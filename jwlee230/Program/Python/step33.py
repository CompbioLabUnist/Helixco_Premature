"""
step33.py: Read & clearify raw TSV for LefSe
"""
import argparse
import pandas
import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("metadata", help="Metadata TSV file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--c", help="Class used for LefSe", type=str, default="Simple Premature")

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")

    raw_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=0).groupby(by="taxonomy").sum()
    print(raw_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1])
    print(sorted(metadata.columns))
    assert args.c in list(metadata.columns)
    print(metadata)

    raw_data.loc["!Subject"] = list(map(lambda x: str(x).replace(" ", "_"), metadata[args.c]))
    raw_data.sort_index(inplace=True)
    print(raw_data)

    for site in tqdm.tqdm(["Cervix", "Mouth", "Vagina", "Neonate-3day", "Neonate-5day"]):
        selected_IDs = metadata.loc[(metadata["Site"] == site), "#SampleID"]
        tmp_data = raw_data.loc[:, selected_IDs]
        tmp_data.to_csv("{0}.{1}.tsv".format(args.output, site), sep="\t")
