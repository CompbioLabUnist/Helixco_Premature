"""
step08.py: draw ANCOM plot
"""
import argparse
import adjustText
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

    if not args.ancom.endswith(".tsv"):
        raise ValueError("ANCOM file must end with .TSV!!")
    elif not args.data.endswith(".tsv"):
        raise ValueError("Data file must end with .TSV!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")

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

    seaborn.scatterplot(data=input_data, x="clr", y="W", ax=ax, style="Reject null hypothesis", markers={True: "o", False: "X"}, legend="full", s=10 ** 3, hue="Reject null hypothesis", palette={True: "tab:blue", False: "tab:red"}, edgecolor=None)

    texts = list()
    for index, row in input_data.loc[(input_data["Reject null hypothesis"]), :].iterrows():
        texts.append(matplotlib.pyplot.text(row["clr"], row["W"], step00.consistency_taxonomy(index, 1), color="black", fontsize="xx-small"))

    matplotlib.pyplot.title("{0} taxa reject null hypothesis".format(input_data.loc[(input_data["Reject null hypothesis"])].shape[0]))
    matplotlib.pyplot.tight_layout()

    adjustText.adjust_text(texts, arrowprops=dict(arrowstyle="-", color="black", alpha=0.3), lim=step00.small, ax=ax)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
