"""
step73.py: statistical test for PTB of mother
"""
import argparse
import collections
import numpy
import pandas
import scipy.stats
import tqdm

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("metadata", help="Metadata TSV file", type=str)
    parser.add_argument("output", help="Output TEX file", type=str)

    args = parser.parse_args()

    if not args.metadata.endswith(".tsv"):
        raise ValueError("METADATA must end with .TSV!!")
    elif not args.output.endswith(".tex"):
        raise ValueError("Output file must end with .TEX!!")

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)
    print(sorted(metadata.columns))

    metadata = metadata.loc[(metadata["Site"] == "Mouth")]
    print(metadata)

    for column in tqdm.tqdm(sorted(metadata.columns)):
        print(column, collections.Counter(metadata[column]).most_common())

    metadata.drop_duplicates(subset=["Data", "Mother"], keep="first", inplace=True)
    print(metadata)

    PTB_data = metadata.loc[(metadata["Gestational Week"] < 34)]
    normal_data = metadata.loc[(metadata["Gestational Week"] >= 34)]

    raw_output_data = list()
    for column in tqdm.tqdm(["Cholesterol", "DBP", "Glucose", "HDL", "Hb", "Hct", "LDL", "Mother Age", "SBP", "Weight gain"]):
        PTB = PTB_data.loc[(PTB_data[column] != -1), column]
        normal = normal_data.loc[(normal_data[column] != -1), column]

        p = scipy.stats.mannwhitneyu(PTB, normal)[1]
        raw_output_data.append([column, f"{numpy.mean(PTB):.1f}±{numpy.std(PTB):.1f}", f"{numpy.mean(normal):.1f}±{numpy.std(normal):.1f}", f"{p:.3f}", "*" if (p < 0.05) else ""])

    for column in tqdm.tqdm(["Advanced maternal age", "C-section", "Gestational Diabetes", "Hypertension", "Mother Antibiotics", "Mother Steroid", "Obesity", "Overweight or Obesity", "PROM", "Preterm Labor", "Too much weight gain"]):
        try:
            PTB = PTB_data[(PTB_data[column])]
            normal = normal_data[(normal_data[column])]
        except KeyError:
            continue

        p = scipy.stats.fisher_exact([[len(PTB), len(PTB_data) - len(PTB)], [len(normal), len(normal_data) - len(normal)]])[1]
        raw_output_data.append([column, f"{len(PTB)} ({len(PTB) / len(PTB_data) * 100:.1f}%)", f"{len(normal)} ({len(normal) / len(normal_data) * 100:.1f}%)", f"{p:.3f}", "*" if (p < 0.05) else ""])

    output_data = pandas.DataFrame(data=raw_output_data, columns=["Clinical", f"<34 GW (n={len(PTB_data)})", f"≥34 GW (n={len(normal_data)})", "p-value", "Remarks"]).set_index("Clinical")
    print(output_data)

    output_data.to_csv(args.output.replace(".tex", ".tsv"), sep="\t")
    output_data.to_latex(args.output, column_format="lrrrc")
