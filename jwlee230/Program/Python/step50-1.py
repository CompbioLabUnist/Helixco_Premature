"""
step50-1.py: Prepare for DESeq2 input
"""
import argparse
import pandas
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tsv file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tsv!!")
    if not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata must end with .tsv!!")

    input_data = pandas.read_csv(args.input, skiprows=1, sep="\t", index_col=["taxonomy", "#OTU ID"])
    input_data.index = list(map(lambda x: f"{x[0]}:{x[1]}", list(input_data.index)))
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col=0)
    metadata["Simple Premature"] = list(map(lambda x: x.replace("+", "_").replace(" ", ""), metadata["Simple Premature"]))
    print(metadata)

    for site in tqdm.tqdm(step00.selected_long_sites):
        selected_metadata = metadata.loc[(metadata["Site"] == site)].sort_values("Simple Premature", ascending=False)
        selected_metadata[["Simple Premature"]].to_csv(f"{args.output}.{site}.coldata", sep="\t")
        input_data.loc[:, list(selected_metadata.index)].to_csv(f"{args.output}.{site}.tsv", sep="\t")
