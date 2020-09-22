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

    parser.add_argument("train", type=str, nargs=1, help="Train TAR.gz file")
    parser.add_argument("normal", type=str, nargs=1, help="Normal TAR.gz file")
    parser.add_argument("meta", type=str, nargs=1, help="Metadata TSV file")
    parser.add_argument("output", type=str, nargs=1, help="Output basename")
    parser.add_argument("--cpu", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if args.cpu < 1:
        raise ValueError("CPU must be greater than zero")
    elif not args.meta[0].endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV")

    train_data = step00.read_pickle(args.train[0])
    normal_data = step00.read_pickle(args.normal[0])
    metadata = pandas.read_csv(args.meta[0], sep="\t", skiprows=[1])

    intersect_columns = sorted(list(set(train_data.columns) & set(normal_data.columns)))
    train_data = train_data[intersect_columns]
    normal_data = normal_data[intersect_columns]

    tsne_data = pandas.DataFrame(sklearn.manifold.TSNE(n_components=2, init="pca", random_state=0, method="exact", n_jobs=args.cpu).fit_transform(pandas.concat([train_data, normal_data], verify_integrity=True)), columns=["TSNE1", "TSNE2"])
    for column in tsne_data.columns:
        tsne_data[column] = sklearn.preprocessing.scale(tsne_data[column])

    train_data["Answer"] = list(metadata["premature"])
    normal_data["Answer"] = "Normal"

    classifier = sklearn.ensemble.RandomForestClassifier(criterion="entropy", max_features=None, n_jobs=args.cpu, random_state=0, bootstrap=False)
    classifier.fit(train_data[intersect_columns], train_data["Answer"])
    feature_importances = classifier.feature_importances_
    best_features = list(map(lambda x: x[1], sorted(zip(feature_importances, intersect_columns), reverse=True)))[:10]

    seaborn.set(context="poster", style="whitegrid")
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
    seaborn.distplot(feature_importances, kde=False, rug=True, axlabel="FeatureImportances")
    fig.savefig(args.output[0] + ".feature_importances.png")
    matplotlib.pyplot.close(fig)

    classifier = sklearn.ensemble.RandomForestClassifier(criterion="entropy", max_features=None, n_jobs=args.cpu, random_state=0, bootstrap=False, warm_start=True)
    for train_index, test_index in sklearn.model_selection.KFold(n_splits=5, random_state=0, shuffle=True).split(train_data):
        classifier.fit(train_data.iloc[train_index][best_features], train_data.iloc[train_index]["Answer"])
        print(classifier.score(train_data.iloc[test_index][best_features], train_data.iloc[test_index]["Answer"]))

    print("Final:", classifier.score(normal_data[best_features], normal_data["Answer"]))
