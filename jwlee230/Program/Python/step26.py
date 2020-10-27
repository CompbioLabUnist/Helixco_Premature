"""
step26.py: RandomForest Classifier
"""
import argparse
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
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

    helixco_data = step00.read_pickle(args.train[0])
    normal_data = step00.read_pickle(args.normal[0])
    premature_data = step00.read_pickle(args.premature[0])
    metadata = pandas.read_csv(args.meta[0], sep="\t", skiprows=[1])

    intersect_columns = sorted(set(helixco_data.columns) & set(normal_data.columns) & set(premature_data.columns))
    helixco_data = helixco_data[intersect_columns]
    normal_data = normal_data[intersect_columns]
    premature_data = premature_data[intersect_columns]

    helixco_data["Answer"] = list(metadata["premature"])
    normal_data["Answer"] = "Normal"
    premature_data["Answer"] = "Premature"

    # Get Feature Importances
    classifier = sklearn.ensemble.RandomForestClassifier(max_features=None, n_jobs=args.cpu, random_state=0)
    classifier.fit(helixco_data[intersect_columns], helixco_data["Answer"])
    feature_importances = classifier.feature_importances_
    best_features = list(map(lambda x: x[1], sorted(list(filter(lambda x: x[0] > 0, zip(feature_importances, intersect_columns))), reverse=True)))

    # Draw Feature Importances
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
    seaborn.distplot(list(filter(lambda x: x > 0, feature_importances)), hist=True, kde=False, rug=True, ax=ax)
    matplotlib.pyplot.title("Feature Importances by Feature Counts")
    matplotlib.pyplot.xlabel("Feature Importances")
    matplotlib.pyplot.ylabel("Counts")
    matplotlib.pyplot.grid(True)
    tar_files.append("importances.png")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    # Select best features
    flag = True
    while flag:
        print(len(best_features), "features!!")

        flag = False
        classifier.fit(helixco_data[best_features], helixco_data["Answer"])
        feature_importances = classifier.feature_importances_
        best_features = list(map(lambda x: x[1], sorted(list(filter(lambda x: x[0] > 0, zip(feature_importances, best_features))), reverse=True)))

        if list(filter(lambda x: x <= 0, feature_importances)):
            flag = True

    # Run K-fold
    k_fold = sklearn.model_selection.StratifiedKFold(n_splits=10)
    test_scores = list()
    for i in range(1, len(best_features) + 1):
        print("With", i, "/", len(best_features), "features!!")
        used_columns = intersect_columns[:i]
        for j, (train_index, test_index) in enumerate(k_fold.split(helixco_data[used_columns], helixco_data["Answer"])):
            x_train, x_test = helixco_data.iloc[train_index][used_columns], helixco_data.iloc[test_index][used_columns]
            y_train, y_test = helixco_data.iloc[train_index]["Answer"], helixco_data.iloc[test_index]["Answer"]

            classifier.fit(x_train, y_train)
            test_scores.append(("Helixco", i, classifier.score(x_test, y_test)))
            test_scores.append(("EBI", i, classifier.score(normal_data[used_columns], normal_data["Answer"])))
            test_scores.append(("HMP", i, classifier.score(premature_data[used_columns], premature_data["Answer"])))

    # Draw K-fold
    score_data = pandas.DataFrame.from_records(test_scores, columns=["Database", "Features", "Accuracy"])
    seaborn.set(context="poster", style="whitegrid")
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
    seaborn.lineplot(data=score_data, x="Features", y="Accuracy", hue="Database", style="Database", ax=ax)
    matplotlib.pyplot.grid(True)
    tar_files.append("accuracy.png")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    # Save data
    with tarfile.open(args.output[0], "w") as tar:
        for file_name in tar_files:
            tar.add(file_name, arcname=file_name)
