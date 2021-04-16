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

    helixco_data: pandas.DataFrame = pandas.concat(list(map(step00.read_pickle, args.helixco)), axis="index", ignore_index=True, copy=False).fillna(value=0)
    EBI_data: pandas.DataFrame = step00.read_pickle(args.ebi)
    HMP_data: pandas.DataFrame = step00.read_pickle(args.hmp)

    same_taxos = sorted(set(helixco_data.columns) & set(EBI_data.columns) & set(HMP_data.columns))

    print(helixco_data.shape, EBI_data.shape, HMP_data.shape)
    print(len(set(helixco_data.columns) & set(EBI_data.columns)), len(set(EBI_data.columns) & set(HMP_data.columns)), len(set(HMP_data.columns) & set(helixco_data.columns)))
    print(len(set(helixco_data.columns) & set(EBI_data.columns) & set(HMP_data.columns)))

    helixco_data = helixco_data[same_taxos]
    EBI_data = EBI_data[same_taxos]
    HMP_data = HMP_data[same_taxos]

    helixco_data["DB"] = "Helixco"
    EBI_data["DB"] = "EBI"
    HMP_data["DB"] = "HMP"

    whole_data = pandas.concat([helixco_data, EBI_data, HMP_data], ignore_index=True, verify_integrity=True)

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="barnes_hut", n_jobs=args.cpu).fit_transform(whole_data[same_taxos]), columns=["tSNE1", "tSNE2"])
    for column in tsne_data.columns:
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])
    tsne_data["DB"] = whole_data["DB"]

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=tsne_data, x="tSNE1", y="tSNE2", ax=ax, hue="DB", style="DB", legend="full", alpha=0.3, s=1000)

    fig.savefig(args.output)
    matplotlib.pyplot.close(fig)
