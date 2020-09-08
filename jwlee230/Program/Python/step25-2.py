"""
step25-2.py: draw t-SNE with pseudo-sample
"""
import argparse
import matplotlib
import matplotlib.pyplot
import seaborn
import pandas
import sklearn.manifold
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("real", type=str, nargs=1, help="Real TAR.gz file")
    parser.add_argument("pseudo", type=str, nargs=1, help="Pseudo TAR.gz file")
    parser.add_argument("output", type=str, nargs=1, help="Output PNG file")
    parser.add_argument("--cpu", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if args.cpu < 1:
        raise ValueError("CPU must be greater than zero")
    elif not args.output[0].endswith(".png"):
        raise ValueError("Output file must end with .PNG")

    real_data = step00.read_pickle(args.real[0])
    pseudo_data = step00.read_pickle(args.pseudo[0])

    used_columns = list(real_data.columns)

    real_data["IsReal"] = "Real"
    pseudo_data["IsReal"] = "Pseudo"

    data = pandas.concat([real_data, pseudo_data], ignore_index=True)

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="exact", n_jobs=args.cpu).fit_transform(data[used_columns]), columns=["TSNE1", "TSNE2"])

    for column in tsne_data.columns:
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])

    tsne_data["IsReal"] = data["IsReal"]

    seaborn.set(context="poster", style="whitegrid")
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

    seaborn.scatterplot(data=tsne_data, x="TSNE1", y="TSNE2", ax=ax, hue="IsReal", style="IsReal", legend="full")

    fig.savefig(args.output[0])
    matplotlib.pyplot.close(fig)
