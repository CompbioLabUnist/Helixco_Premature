"""
step22-2.py: draw t-SNE with multiple databases
"""
import argparse
import pandas
import matplotlib
import matplotlib.pyplot
import seaborn
import sklearn.manifold
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("helixco", type=str, nargs="+", help="Helixco TAR.gz file")
    parser.add_argument("ebi", type=str, help="EBI TAR.gz file")
    parser.add_argument("hmp", type=str, help="HMP TAR.gz file")
    parser.add_argument("output", type=str, help="Output PDF file")
    parser.add_argument("--cpu", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    if args.cpu < 1:
        raise ValueError("CPU must be greater than zero!!")

    helixco_data: pandas.DataFrame = pandas.concat(list(map(step00.read_pickle, args.helixco)), axis="index", join="inner", ignore_index=True, copy=False)
    EBI_data: pandas.DataFrame = step00.read_pickle(args.ebi)
    HMP_data: pandas.DataFrame = step00.read_pickle(args.hmp)

    same_taxa = set(helixco_data.columns) & set(EBI_data.columns) & set(HMP_data.columns)
    two_taxa = (set(helixco_data.columns) & set(EBI_data.columns)) - same_taxa, (set(EBI_data.columns) & set(HMP_data.columns)) - same_taxa, (set(HMP_data.columns) & set(helixco_data.columns)) - same_taxa
    one_taxa = set(helixco_data.columns) - set(EBI_data.columns) - set(HMP_data.columns), set(EBI_data.columns) - set(HMP_data.columns) - set(helixco_data.columns), set(HMP_data.columns) - set(helixco_data.columns) - set(EBI_data.columns)

    print(list(map(len, one_taxa)))
    print(list(map(len, two_taxa)))
    print(len(same_taxa))

    helixco_data = helixco_data[sorted(same_taxa)]
    EBI_data = EBI_data[sorted(same_taxa)]
    HMP_data = HMP_data[sorted(same_taxa)]

    helixco_data["DB"] = "Helixco"
    EBI_data["DB"] = "EBI"
    HMP_data["DB"] = "HMP"

    whole_data = pandas.concat([EBI_data, HMP_data, helixco_data], join="inner", ignore_index=True, verify_integrity=True)

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="barnes_hut", n_jobs=args.cpu).fit_transform(whole_data[same_taxa]), columns=["tSNE1", "tSNE2"])
    for column in tsne_data.columns:
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])
    tsne_data["DB"] = whole_data["DB"]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    order = ["Helixco", "EBI", "HMP"]
    seaborn.scatterplot(data=tsne_data, x="tSNE1", y="tSNE2", ax=ax, hue="DB", style="DB", hue_order=order, style_order=order, legend="full", alpha=0.3, s=1000, edgecolor="none")
    matplotlib.pyplot.title("Total {0} taxa and {1} samples".format(whole_data.shape[1], whole_data.shape[0]))

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
