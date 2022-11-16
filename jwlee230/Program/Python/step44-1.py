"""
step44-1.py: make DESeq2 input TSV
"""
import argparse
import pandas
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("output", help="Output file basename", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=["#OTU ID"]).groupby("taxonomy").sum().T
    input_data = input_data.loc[:, list(filter(step00.filtering_taxonomy, list(input_data.columns)))].T
    input_data.index.name = "gene_name"
    print(input_data)

    for key, value in tqdm.tqdm(step00.selected_sites_dict.items()):
        tmp_data = input_data.loc[:, list(filter(lambda x: x.endswith(key), list(input_data.columns)))]
        tmp_data.to_csv(f"{args.output}.{value}.tsv", sep="\t")
