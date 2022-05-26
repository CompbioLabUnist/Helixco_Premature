"""
step48.py: Shared taxonomy proportion with neonates and mother
"""
import argparse
import itertools
import multiprocessing
import typing
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import statannotations.Annotator
import step00

input_data = pandas.DataFrame()


def corresponding(x: str) -> typing.List[typing.Tuple[str, str, str, str]]:
    xs = x.split("-")
    answer = []

    if (i := "{0}-{1}-Mother-C".format(xs[0], xs[1])) in list(input_data.columns):
        answer.append((x, i, xs[3], "C"))
    if (i := "{0}-{1}-Mother-M".format(xs[0], xs[1])) in list(input_data.columns):
        answer.append((x, i, xs[3], "M"))
    if (i := "{0}-{1}-Mother-V".format(xs[0], xs[1])) in list(input_data.columns):
        answer.append((x, i, xs[3], "V"))

    return answer


def calculate(x: str, y: str) -> float:
    baby = set(input_data.loc[(input_data[x] != 0)].index)
    mother = set(input_data.loc[(input_data[y] != 0)].index)
    return len(baby & mother) / len(baby | mother)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .tar.gz!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("Number of CPUs must be positive!!")

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col="#OTU ID")
    del input_data["taxonomy"]
    print(input_data)

    baby_index = sorted(filter(lambda x: (x.split("-")[-1] in {"B1", "B3", "B5"}), list(input_data.columns)))
    print(baby_index)

    with multiprocessing.Pool(args.cpus) as pool:
        output_data = pandas.DataFrame(columns=["Baby", "Mother", "BabyType", "MotherType"], data=itertools.chain.from_iterable(pool.map(corresponding, baby_index)))
        output_data["Shared"] = pool.starmap(calculate, output_data[["Baby", "Mother"]].itertuples(index=False, name=None))
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(30, 24))

    seaborn.boxplot(data=output_data, x="BabyType", y="Shared", order=["B1", "B3", "B5"], hue="MotherType", hue_order=["C", "M", "V"], ax=ax)
    statannotations.Annotator.Annotator(ax, [((b1, m), (b2, m)) for m in ["C", "M", "V"] for b1, b2 in list(itertools.combinations(["B1", "B3", "B5"], 2))], data=output_data, x="BabyType", y="Shared", order=["B1", "B3", "B5"], hue="MotherType", hue_order=["C", "M", "V"]).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()

    matplotlib.pyplot.xlabel("Sample Type")
    matplotlib.pyplot.ylabel("Shared Taxonomy Proportion")
    matplotlib.pyplot.tight_layout()

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
