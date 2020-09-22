"""
step26.py: RandomForest Classifier
"""
import argparse
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import scipy.stats
import seaborn
import sklearn.ensemble
import sklearn.manifold
import sklearn.model_selection
import sklearn.tree
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
        raise ValueError("CPU must be greater than zero!!")
    elif not args.meta[0].endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output[0].endswith(".tar"):
        raise ValueError("Output file must end with .tar!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 30})

    tar_files: typing.List[str] = list()

    train_data = step00.read_pickle(args.train[0])
    normal_data = step00.read_pickle(args.normal[0])
    metadata = pandas.read_csv(args.meta[0], sep="\t", skiprows=[1])

    intersect_columns = sorted(list(set(train_data.columns) & set(normal_data.columns)))
    train_data = train_data[intersect_columns]
    normal_data = normal_data[intersect_columns]

    train_data["Answer"] = list(metadata["premature"])
    normal_data["Answer"] = "Normal"

    # Get Feature Importances
    classifier = sklearn.ensemble.RandomForestClassifier(criterion="entropy", max_features=None, n_jobs=args.cpu, random_state=0, bootstrap=False)
    classifier.fit(train_data[intersect_columns], train_data["Answer"])
    feature_importances = classifier.feature_importances_
    best_features = list(map(lambda x: x[1], sorted(list(filter(lambda x: x[0] > 0, zip(feature_importances, intersect_columns))), reverse=True)))

    # Draw Feature Importances
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
    seaborn.distplot(list(filter(lambda x: x > 0, feature_importances)), hist=True, kde=False, rug=True, ax=ax)
    matplotlib.pyplot.title("Feature Importances")
    matplotlib.pyplot.xlabel("Feature")
    matplotlib.pyplot.ylabel("Importances")
    tar_files.append("importances.png")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(train_data[best_features], train_data["Answer"], test_size=0.1, random_state=0, shuffle=True, stratify=train_data["Answer"])

    # Select best features
    scores = [0 for _ in best_features]
    best_num, best_score = 0, 0
    for i in range(1, len(scores) + 1):
        classifier.fit(x_train[best_features[:i]], y_train)
        scores[i - 1] = tmp = classifier.score(x_test[best_features[:i]], y_test)
        if tmp > best_score:
            best_score = tmp
            best_num = i

    best_features = best_features[:best_num]
    x_train = x_train[best_features]
    x_test = x_test[best_features]
    for f in best_features:
        print(f)

    # Draw scores
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
    matplotlib.pyplot.plot(range(1, len(scores) + 1), scores, "o-")
    matplotlib.pyplot.xlabel("Number of Features")
    matplotlib.pyplot.ylabel("Accuracy")
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.title("Best Features: %d Features at %.2f" % (best_num, best_score))
    tar_files.append("scores.png")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    # Draw Random Forest
    fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
    tree = classifier.fit(x_train, y_train)
    sklearn.tree.plot_tree(tree.estimators_[0], ax=ax, feature_names=[list(filter(lambda x: not x.endswith("__"), x))[-1].strip() for x in list(map(lambda x: x.split(";"), best_features))], class_names=["Normal", "Premature"], filled=True)
    matplotlib.pyplot.title("Accuracy: %.2f" % best_score)
    tar_files.append("tree.png")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    # Draw violin plot
    for i, feature in enumerate(best_features):
        data = train_data[["Answer"] + [feature]]
        seaborn.set(context="poster", style="whitegrid")
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
        seaborn.violinplot(data=data, x="Answer", y=feature, ax=ax, scale="count", inner="stick")
        matplotlib.pyplot.xlabel("Normal/Premature")
        matplotlib.pyplot.ylabel(list(filter(lambda x: not x.endswith("__"), feature.split(";")))[-1])
        matplotlib.pyplot.title("P-value: %.2f" % scipy.stats.ttest_ind(data.loc[(data["Answer"] == "Normal")][feature], data.loc[(data["Answer"] == "Premature")][feature], equal_var=False)[1])
        tar_files.append("feature_%d.png" % i)
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output[0], "w") as tar:
        for f in tar_files:
            tar.add(f, arcname=f)
