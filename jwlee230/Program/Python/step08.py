"""
step08.py: draw ANCOM plot
"""
import argparse
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("ancom", help="ANCOM ancom.tsv file", type=str)
    parser.add_argument("data", help="ANCOM data.tsv file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)

    args = parser.parse_args()

    ancom_data = pandas.read_csv(args.ancom, sep="\t", names=["id", "W", "Reject null hypothesis"], usecols=["id", "Reject null hypothesis"], header=0, index_col="id")
    print(ancom_data)

    input_data = pandas.read_csv(args.data, sep="\t", index_col="id")
    print(input_data)

    input_data = pandas.concat([input_data, ancom_data], axis="columns", verify_integrity=True)
    print(input_data)

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

    seaborn.scatterplot(data=input_data, x="clr", y="W", ax=ax, style="Reject null hypothesis", markers={True: "o", False: "X"}, legend="full", s=1000, hue="Reject null hypothesis", palette={True: "tab:blue", False: "tab:red"})

    matplotlib.pyplot.title("{0} taxa reject null hypothesis".format(input_data.loc[(input_data["Reject null hypothesis"])].shape[0]))

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
