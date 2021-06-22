"""
step22-3.py: draw t-SNE plot
"""
import argparse
import itertools
import tarfile
import pandas
import matplotlib
import matplotlib.pyplot
import seaborn
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Input TAR.gz file")
    parser.add_argument("output", type=str, help="Output TAR file")

    args = parser.parse_args()

    if not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .tar!!")

    tar_files = list()

    input_data: pandas.DataFrame = step00.read_pickle(args.input)
    print(input_data)

    columns = set(input_data.columns)
    columns -= {"tSNE1", "tSNE2"}
    print(sorted(columns))

    numeric_columns = step00.numeric_columns & columns
    categorical_columns = columns - step00.numeric_columns
    print(numeric_columns)
    print(categorical_columns)

    for c in numeric_columns:
        print("Numeric:", c)
        matplotlib.use("Agg")
        matplotlib.rcParams.update(step00.matplotlib_parameters)
        seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
        seaborn.scatterplot(data=input_data, x="tSNE1", y="tSNE2", ax=ax, hue=c, size=c, legend="brief")

        tar_files.append(c.replace(" ", "_") + ".pdf")
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    for c in categorical_columns:
        print("Categorical:", c)
        matplotlib.use("Agg")
        matplotlib.rcParams.update(step00.matplotlib_parameters)
        seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
        seaborn.scatterplot(data=input_data, x="tSNE1", y="tSNE2", ax=ax, hue=c, style=c, legend="full", markers={x: y for x, y in zip(sorted(set(input_data[c])), itertools.cycle(step00.markers))}, s=300, style_order=sorted(set(input_data[c])), hue_order=sorted(set(input_data[c])))

        tar_files.append(c.replace(" ", "_") + ".pdf")
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for file_name in tar_files:
            tar.add(file_name, arcname=file_name)
