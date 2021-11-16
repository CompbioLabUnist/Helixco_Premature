"""
step39-1.py: Beta-diversity heatmap plot for Macrogen data
"""
import argparse
import multiprocessing
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00

input_data = pandas.DataFrame()


def corresponding(x: str) -> str:
    xs = x.split("-")
    return "{0}-{1}-Mother-V".format(xs[0], xs[1])


def calculate(x: str) -> float:
    baby = set(input_data.loc[(input_data[x] != 0)].index)
    mother = set(input_data.loc[(input_data[corresponding(x)] != 0)].index)
    return len(baby & mother) / len(baby)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of CPUs", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tar.gz"):
        raise ValueError("Input file must end with .tar.gz!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("Number of CPUs must be positive!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = step00.read_pickle(args.input)
    input_data.index = list(map(step00.simplified_taxonomy, list(input_data.index)))
    input_data = input_data.iloc[:, 1:]
    print(input_data)

    baby_index = list(filter(lambda x: (x.split("-")[-1] in {"B1", "B3", "B5"}) and (corresponding(x) in set(input_data.columns)), list(input_data.columns)))
    print(baby_index)

    output_data = pandas.DataFrame(index=baby_index)
    output_data["Type"] = list(map(lambda x: x.split("-")[-1], list(output_data.index)))
    with multiprocessing.Pool(args.cpus) as pool:
        output_data["Shared"] = pool.map(calculate, list(output_data.index))
    print(output_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.violinplot(data=output_data, x="Type", y="Shared", order=["B1", "B3", "B5"], inner="box", ax=ax)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
