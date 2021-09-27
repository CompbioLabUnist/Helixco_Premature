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
import tqdm
import step00

data = pandas.DataFrame()


def read(file_name: str) -> pandas.DataFrame:
    return pandas.read_csv(file_name, sep="\t", index_col=0)


def draw(alpha: str, disease: str, site: str) -> str:
    print(alpha, disease, site)

    drawing_data = data.loc[(data["Site"] == site)]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
    seaborn.violinplot(data=drawing_data, x=disease, y=alpha, order=sorted(set(data[disease])), inner="box", ax=ax)

    try:
        statannot.add_stat_annotation(ax, data=drawing_data, x=disease, y=alpha, order=sorted(set(data[disease])), test="Mann-Whitney", box_pairs=itertools.combinations(sorted(set(data[disease])), 2), text_format="simple", loc="outside", verbose=0, fontsize=step00.matplotlib_parameters["font.size"])
    except ValueError:
        pass

    matplotlib.pyplot.ylabel(alpha.replace("_", " "))
    matplotlib.pyplot.tight_layout()

    file_name = "{2}+{1}+{0}.pdf".format(alpha, disease.replace(" ", "_"), site)
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


def draw_all(alpha: str, disease: str) -> str:
    print(alpha, disease)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
    seaborn.violinplot(data=data, x=disease, y=alpha, order=sorted(set(data[disease])), inner="box", ax=ax)

    try:
        statannot.add_stat_annotation(ax, data=data, x=disease, y=alpha, order=sorted(set(data[disease])), test="Mann-Whitney", box_pairs=itertools.combinations(sorted(set(data[disease])), 2), text_format="simple", loc="outside", verbose=0, fontsize=step00.matplotlib_parameters["font.size"])
    except ValueError:
        pass

    matplotlib.pyplot.ylabel(alpha.replace("_", " "))
    matplotlib.pyplot.tight_layout()

    file_name = "All+{1}+{0}.pdf".format(alpha, disease.replace(" ", "_"))
    fig.savefig(file_name)
    matplotlib.pyplot.close(fig)
    return file_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TSV file", nargs="+")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output TAR file")
    parser.add_argument("--cpus", type=int, default=1, help="CPU to use")

    data_group = parser.add_mutually_exclusive_group()
    data_group.add_argument("--first", help="Select First data", action="store_true", default=False)
    data_group.add_argument("--second", help="Select Second+Third data", action="store_true", default=False)

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
    alphas = set(raw_data.columns)

    if args.first:
        raw_data = raw_data.loc[list(filter(lambda x: x.startswith("First"), list(raw_data.index))), :]
    elif args.second:
        raw_data = raw_data.loc[list(filter(lambda x: x.startswith("Second") or x.startswith("Third"), list(raw_data.index))), :]

    print(raw_data)
    print(sorted(alphas))

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1], dtype=str).dropna(axis="columns", how="all").set_index(keys=["#SampleID"], verify_integrity=True)
    metadata = metadata.loc[list(raw_data.index), sorted(set(metadata.columns) - step00.numeric_columns)].replace(to_replace=-1, value=None)
    diseases = set(metadata.columns) - step00.numeric_columns - {"Mother", "Neonate"}
    print(metadata)
    print(sorted(diseases))

    data = pandas.concat(objs=[raw_data, metadata], axis="columns", verify_integrity=True)
    sites = set(data["Site"])
    print(data)
    print(sorted(sites))

    with multiprocessing.Pool(args.cpus) as pool:
        files = pool.starmap(draw, itertools.product(alphas, diseases, sites))
        files += pool.starmap(draw_all, itertools.product(alphas, diseases))

    with tarfile.open(args.output, "w") as tar:
        for f in tqdm.tqdm(files):
            tar.add(f, arcname=f)
