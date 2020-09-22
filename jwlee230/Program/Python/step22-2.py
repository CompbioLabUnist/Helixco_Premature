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

    parser.add_argument("helixco", type=str, nargs=1, help="Helixco TAR.gz file")
    parser.add_argument("ebi", type=str, nargs=1, help="EBI TAR.gz file")
    parser.add_argument("output", type=str, nargs=1, help="Output PNG file")
    parser.add_argument("--cpu", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if args.cpu < 1:
        raise ValueError("CPU must be greater than zero!!")

    helixco_data: pandas.DataFrame = step00.read_pickle(args.helixco[0])
    EBI_data: pandas.DataFrame = step00.read_pickle(args.ebi[0])

    same_taxos = sorted(set(helixco_data.columns) & set(EBI_data.columns))

    helixco_data = helixco_data[same_taxos]
    EBI_data = EBI_data[same_taxos]

    helixco_data["DB"] = "Helixco"
    EBI_data["DB"] = "EBI"

    whole_data = pandas.concat([helixco_data, EBI_data], ignore_index=True, verify_integrity=True)

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="exact", n_jobs=args.cpu).fit_transform(whole_data[same_taxos]), columns=["TSNE1", "TSNE2"])
    for column in tsne_data.columns:
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])
    tsne_data["DB"] = whole_data["DB"]

    seaborn.set(context="poster", style="whitegrid")
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=tsne_data, x="TSNE1", y="TSNE2", ax=ax, hue="DB", style="DB", legend="full")

    fig.savefig(args.output[0])
    matplotlib.pyplot.close(fig)
