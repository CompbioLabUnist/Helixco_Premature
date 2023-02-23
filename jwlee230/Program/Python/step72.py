"""
step72.py: PICRUSt2 results to TSV
"""
import argparse
import itertools
import multiprocessing
import numpy
import pandas
import scipy.stats
import tqdm
import step00

input_data = pandas.DataFrame()


def test(pathway, a_list, b_list):
    p = scipy.stats.ttest_ind(input_data.loc[pathway, a_list], input_data.loc[pathway, b_list])[1]
    if p > 0.05:
        return ""

    if numpy.mean(input_data.loc[pathway, a_list]) < numpy.mean(input_data.loc[pathway, b_list]):
        return f"Up ({p:.3f})"
    else:
        return f"Down ({p:.3f})"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV.gz file", type=str)
    parser.add_argument("metadata", help="Metadata TSV file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs to use", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tsv.gz"):
        raise ValueError("INPUT must end with .TSV.gz!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("METADATA must end with .TSV!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output file must end with .TSV!!")
    elif args.cpus < 1:
        raise ValueError("Number of CPUs must be positive!!")

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0)
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    output_data = pandas.DataFrame(index=input_data.index)
    output_data["description"] = input_data["description"]
    print(output_data)

    pathways = list(input_data.index)
    for site, (PTB1, PTB2) in tqdm.tqdm(list(itertools.product(step00.selected_long_sites, itertools.combinations(("PTB", "Normal"), r=2)))):
        a_list = list(metadata.loc[(metadata["Site"] == site) & (metadata["Premature"] == PTB1)].index)
        b_list = list(metadata.loc[(metadata["Site"] == site) & (metadata["Premature"] == PTB2)].index)

        with multiprocessing.Pool(args.cpus) as pool:
            output_data[f"{site}: {PTB1}-{PTB2}"] = pool.starmap(test, [(pathway, a_list, b_list) for pathway in pathways])
    print(output_data)
    output_data.to_csv(args.output, sep="\t")
