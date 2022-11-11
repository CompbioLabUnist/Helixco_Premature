"""
step11.py: making t-SNE
"""
import argparse
import pandas
import sklearn.manifold
import sklearn.preprocessing
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output TSV file")
    parser.add_argument("--cpus", type=int, default=1, help="CPUs to use")

    args = parser.parse_args()

    if args.cpus < 1:
        raise ValueError("CPUS must be greater than zero")
    elif not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file muste end with .TSV!!")

    raw_data = pandas.read_csv(args.input, sep="\t", skiprows=1)
    raw_data.set_index(inplace=True, keys=["taxonomy", "#OTU ID"], verify_integrity=True)
    raw_data = raw_data.T
    print(raw_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[sorted(set(raw_data.index) & set(metadata.index)), :].replace(to_replace=-1, value=None)
    print(metadata)

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=42, method="exact", n_jobs=args.cpus, perplexity=50, n_iter=step00.big, verbose=1).fit_transform(raw_data), columns=["tSNE1", "tSNE2"])

    for column in tqdm.tqdm(list(tsne_data.columns)):
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])
        metadata[column] = list(tsne_data[column])
    print(metadata)

    metadata.to_csv(args.output, sep="\t")
