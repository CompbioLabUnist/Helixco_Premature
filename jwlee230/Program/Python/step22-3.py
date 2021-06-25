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
    print(sorted(numeric_columns))
    print(sorted(categorical_columns))

    for site in set(input_data["Site"]):
        data = input_data.loc[(input_data["Site"] == site)]

        for c in numeric_columns:
            print("Numeric:", c)
            matplotlib.use("Agg")
            matplotlib.rcParams.update(step00.matplotlib_parameters)
            seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
            seaborn.scatterplot(data=data, x="tSNE1", y="tSNE2", ax=ax, hue="Detail Premature", size=c, legend="brief", hue_order=step00.detailed_PTB)

            tar_files.append("{0}+{1}.pdf".format(site, c.replace(" ", "_")))
            fig.savefig(tar_files[-1])
            matplotlib.pyplot.close(fig)

        for c in categorical_columns:
            print("Categorical:", c)
            matplotlib.use("Agg")
            matplotlib.rcParams.update(step00.matplotlib_parameters)
            seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
            seaborn.scatterplot(data=data, x="tSNE1", y="tSNE2", ax=ax, hue="Detail Premature", style=c, legend="brief", markers={x: y for x, y in zip(sorted(set(input_data[c])), itertools.cycle(step00.markers))}, s=40 ** 2, style_order=sorted(set(input_data[c])), hue_order=step00.detailed_PTB)

            tar_files.append("{0}+{1}.pdf".format(site, c.replace(" ", "_")))
            fig.savefig(tar_files[-1])
            matplotlib.pyplot.close(fig)

    for c in numeric_columns:
        print("Numeric:", c)
        matplotlib.use("Agg")
        matplotlib.rcParams.update(step00.matplotlib_parameters)
        seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
        seaborn.scatterplot(data=input_data, x="tSNE1", y="tSNE2", ax=ax, hue="Detail Premature", size=c, legend="brief", hue_order=step00.detailed_PTB)

        tar_files.append("{0}+{1}.pdf".format("Whole", c.replace(" ", "_")))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    for c in categorical_columns:
        print("Categorical:", c)
        matplotlib.use("Agg")
        matplotlib.rcParams.update(step00.matplotlib_parameters)
        seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
        seaborn.scatterplot(data=input_data, x="tSNE1", y="tSNE2", ax=ax, hue="Detail Premature", style=c, legend="brief", markers={x: y for x, y in zip(sorted(set(input_data[c])), itertools.cycle(step00.markers))}, s=40 ** 2, style_order=sorted(set(input_data[c])), hue_order=step00.detailed_PTB)

        tar_files.append("{0}+{1}.pdf".format("Whole", c.replace(" ", "_")))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output, "w") as tar:
        for file_name in tar_files:
            tar.add(file_name, arcname=file_name)
