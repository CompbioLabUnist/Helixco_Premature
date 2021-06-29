"""
step30.py: draw alpha-diversity violin plots
"""
import argparse
import itertools
import multiprocessing
import tarfile
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import statannot
import step00

data = pandas.DataFrame()


def read(file_name: str) -> pandas.DataFrame:
    return pandas.read_csv(file_name, sep="\t", index_col=0)


def draw(alpha: str, disease: str) -> str:
    print(alpha, disease)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(36, 36))
    seaborn.violinplot(data=data, x=disease, y=alpha, hue="Premature", hue_order=sorted(set(data["Premature"])), order=sorted(set(data[disease])), inner="box", ax=ax)

    statannot.add_stat_annotation(ax, data=data, x=disease, y=alpha, hue="Premature", order=sorted(set(data[disease])), test="t-test_ind", box_pairs=itertools.combinations(itertools.product(sorted(set(data[disease])), sorted(set(data["Premature"]))), 2), text_format="star", loc="outside", verbose=2, fontsize=step00.matplotlib_parameters["font.size"])

    matplotlib.pyplot.ylabel(alpha.replace("_", " "))
    matplotlib.pyplot.tight_layout()

    file_name = "{0}+{1}.pdf".format(alpha, disease.replace(" ", "_"))
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file", nargs="+")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output TAR file")
    parser.add_argument("--cpus", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".tsv"), args.input)):
        raise ValueError("Input must end with .TSV!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .TAR!!")
    elif args.cpus < 1:
        raise ValueError("CPUS must be greater than zero!!")

    raw_data = pandas.concat(objs=list(map(read, args.input)), axis="columns", verify_integrity=True)
    alphas = list(raw_data.columns)
    print(raw_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[list(raw_data.index), sorted(set(metadata.columns) - step00.numeric_columns)].replace(to_replace=-1, value=None)
    diseases = list(metadata.columns)
    print(metadata)

    data = pandas.concat(objs=[raw_data, metadata], axis="columns", verify_integrity=True)
    data = data.loc[(data["Site"].isin(step00.selected_sites))]
    print(data)

    with multiprocessing.Pool(args.cpus) as pool:
        files = pool.starmap(draw, itertools.product(alphas, diseases))

    with tarfile.open(args.output, "w") as tar:
        for f in files:
            print("Compressing:", f)
            tar.add(f, arcname=f)
