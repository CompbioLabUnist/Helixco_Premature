"""
step44-2.py: make DESeq2 input coldata
"""
import argparse
import pandas
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("metadata", help="Metadata file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)
    parser.add_argument("--compare", help="Compare clinical data", nargs=3, default=["Simple Premature", "Late PTB+Normal", "Early PTB"])

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata must end with .tsv!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=["#OTU ID"]).groupby("taxonomy").sum().T
    input_data = input_data.loc[:, list(filter(step00.filtering_taxonomy, list(input_data.columns)))]
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], index_col=0)
    print(metadata)
    print(sorted(metadata.columns))

    metadata = metadata.loc[list(filter(lambda x: x in (set(input_data.index) & set(metadata.index)), list(metadata.index)))]
    print(metadata)

    for site in tqdm.tqdm(step00.selected_long_sites):
        tmp_data = metadata.loc[(metadata["Site"] == site)]

        index = list(filter(lambda x: (x in set(tmp_data.loc[(tmp_data[args.compare[0]] == args.compare[1])].index)), list(input_data.index))) + list(filter(lambda x: (x in set(tmp_data.loc[(tmp_data[args.compare[0]] == args.compare[2])].index)), list(input_data.index)))

        output_data = input_data.loc[index, :]
        output_data.T.to_csv(f"{args.output}.{site}.tsv", sep="\t")

        output_data = pandas.DataFrame(index=index)
        output_data["condition"] = list(map(lambda x: metadata.loc[x, args.compare[0]], list(output_data.index)))
        output_data.index.name = "ID"
        output_data.to_csv(f"{args.output}.{site}.coldata", sep="\t")
