"""
step45-1.py: get alpha-diversity
"""
import argparse
import io
import pandas
import skbio
import skbio.diversity
import tqdm
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("tree", type=str, help="Tree NWk file")
    parser.add_argument("output", type=str, help="Output TAR file")

    args = parser.parse_args()

    if not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.input.endswith(".tsv"):
        raise ValueError("Input must end with .TSV!!")
    elif not args.tree.endswith(".nwk"):
        raise ValueError("Tree must end with .NWK!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

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

    output_data = pandas.DataFrame(index=metadata.index)
    alphas = ["faith_pd", "observed_otus", "pielou_e", "shannon"]
    for alpha in tqdm.tqdm(alphas):
        if alpha.endswith("_ci"):
            continue

        try:
            if alpha == "faith_pd":
                output_data[alpha] = skbio.diversity.alpha_diversity(alpha, input_data, ids=list(input_data.index), otu_ids=list(input_data.columns), tree=tree)
            else:
                output_data[alpha] = skbio.diversity.alpha_diversity(alpha, input_data, ids=list(input_data.index))
        except (TypeError, KeyError, ValueError) as e:
            print(alpha, e)

    print(output_data)
    output_data.to_csv(args.output, sep="\t")
