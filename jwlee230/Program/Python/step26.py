"""
step26.py: RandomForest Classifier
"""
import argparse
import matplotlib
import matplotlib.pyplot
import seaborn
import sklearn.ensemble
import sklearn.manifold
import sklearn.model_selection
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("real", type=str, nargs=1, help="Real TAR.gz file")
    parser.add_argument("meta", type=str, nargs=1, help="Metadata TSV file")
    parser.add_argument("output", type=str, nargs=1, help="Output basename")
    parser.add_argument("--cpu", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if args.cpu < 1:
        raise ValueError("CPU must be greater than zero")
    elif not args.meta[0].endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV")

    real_data = step00.read_pickle(args.real[0])
    metadata = pandas.read_csv(args.meta[0], sep="\t", skiprows=[1])
    answer_column = "premature"

    metadata.set_index("#SampleID", inplace=True, verify_integrity=True)

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="exact", n_jobs=args.cpu).fit_transform(real_data), columns=["TSNE1", "TSNE2"])
    for column in tsne_data.columns:
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])
    tsne_data.index = real_data.index

    classifier = sklearn.ensemble.RandomForestClassifier(criterion="entropy", max_features=None, n_jobs=args.cpu, random_state=0)
    classifier.fit(real_data, metadata[answer_column])
    feature_importances = classifier.feature_importances_
    best_features = list(map(lambda x: x[1], sorted(zip(feature_importances, real_data.columns), reverse=True)))[:10]

    seaborn.set(context="poster", style="whitegrid")
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
    seaborn.distplot(feature_importances, kde=False, rug=True, axlabel="FeatureImportances")
    fig.savefig(args.output[0] + ".feature_importances.png")
    matplotlib.pyplot.close(fig)

    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(real_data[best_features], metadata[answer_column], test_size=0.1, random_state=0, shuffle=True, stratify=metadata[answer_column])
    classifier.fit(x_train, y_train)
    prediction = classifier.predict(x_test)

    print(classifier.score(x_test, y_test), list(prediction))
