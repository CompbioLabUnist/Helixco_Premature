"""
step46-1.py: beta-diversity
"""
import argparse
import io
import pandas
import skbio
import skbio.diversity
import skbio.stats
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("tree", type=str, help="Tree NWk file")
    parser.add_argument("output", type=str, help="Output file basename")

    args = parser.parse_args()

    if not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.tree.endswith(".nwk"):
        raise ValueError("Tree must end with .NWK!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID")
    del input_data["taxonomy"]
    input_data = input_data.T
    print(input_data)

    with open(args.tree, "r") as f:
        tree = skbio.TreeNode.read(io.StringIO(f.readline()))

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[sorted(set(input_data.index) & set(metadata.index)), :].replace(to_replace=-1, value=None)
    metadata_columns = sorted(set(metadata.columns) - step00.numeric_columns - {"Mother", "Too much weight gain", "Site"})
    print(metadata)
    print(metadata_columns)

    for beta in skbio.diversity.get_beta_diversity_metrics():
        beta_diversity = skbio.diversity.beta_diversity(beta, input_data, ids=list(input_data.index), otu_ids=list(input_data.columns), tree=tree).to_data_frame()
        print(beta_diversity)
        beta_diversity.to_csv(f"{args.output}.{beta}.tsv", sep="\t")
