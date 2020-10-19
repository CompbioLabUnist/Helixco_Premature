"""
step26.py: Decision Tree Classifier
"""
import argparse
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import scipy.stats
import seaborn
import sklearn.manifold
import sklearn.model_selection
import sklearn.tree
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("train", type=str, nargs=1, help="Train TAR.gz file")
    parser.add_argument("normal", type=str, nargs=1, help="Normal TAR.gz file")
    parser.add_argument("premature", type=str, nargs=1, help="premature TAR.gz file")
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
    premature_data = step00.read_pickle(args.premature[0])
    metadata = pandas.read_csv(args.meta[0], sep="\t", skiprows=[1])

    intersect_columns = sorted(set(train_data.columns) & set(normal_data.columns) & set(premature_data.columns))
    train_data = train_data[intersect_columns]
    normal_data = normal_data[intersect_columns]
    premature_data = premature_data[intersect_columns]

    train_data["Answer"] = list(metadata["premature"])
    normal_data["Answer"] = "Normal"
    premature_data["Answer"] = "Premature"

    # Get Feature Importances
    classifier = sklearn.tree.DecisionTreeClassifier(max_features=None, random_state=0)
    classifier.fit(train_data[intersect_columns], train_data["Answer"])
    feature_importances = classifier.feature_importances_
    best_features = list(map(lambda x: x[1], sorted(list(filter(lambda x: x[0] > 0, zip(feature_importances, intersect_columns))), reverse=True)))

    # Draw Feature Importances
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
    seaborn.distplot(list(filter(lambda x: x > 0, feature_importances)), hist=True, kde=False, rug=True, ax=ax)
    matplotlib.pyplot.title("Feature Importances")
    matplotlib.pyplot.xlabel("Feature Importances")
    matplotlib.pyplot.ylabel("Counts")
    tar_files.append("importances.png")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    validation_data = pandas.concat([normal_data, premature_data], verify_integrity=True)

    # Select best features
    scores = [0 for _ in best_features]
    best_num, best_score = 0, 0
    for i in range(1, len(scores) + 1):
        classifier.fit(train_data[best_features[:i]], train_data["Answer"])
        scores[i - 1] = tmp = classifier.score(validation_data[best_features[:i]], validation_data["Answer"])
        if tmp > best_score:
            best_score = tmp
            best_num = i

    best_features = best_features[:best_num]
    tar_files.append("features.txt")
    with open(tar_files[-1], "w") as f:
        for feature in best_features:
            f.write(feature)
            f.write("\n")

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

    # Draw violin plot
    for i, feature in enumerate(best_features):
        data = pandas.concat([train_data, normal_data, premature_data], verify_integrity=True)[["Answer"] + [feature]]
        seaborn.set(context="poster", style="whitegrid")
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
        seaborn.violinplot(data=data, x="Answer", y=feature, ax=ax)
        matplotlib.pyplot.xlabel("Normal/Premature")
        matplotlib.pyplot.ylabel(list(filter(lambda x: not x.endswith("__"), feature.split(";")))[-1])
        matplotlib.pyplot.title("P-value: %.4f" % scipy.stats.ttest_ind(data.loc[(data["Answer"] == "Normal")][feature], data.loc[(data["Answer"] == "Premature")][feature], equal_var=False)[1])
        tar_files.append("feature_%d.png" % i)
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    with tarfile.open(args.output[0], "w") as tar:
        for file_name in tar_files:
            tar.add(file_name, arcname=file_name)