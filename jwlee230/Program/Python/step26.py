"""
step26.py: RandomForest Classifier
"""
import argparse
import itertools
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import seaborn
import sklearn.ensemble
import sklearn.manifold
import sklearn.model_selection
import sklearn.tree
import statannot
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Train TAR.gz file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output basename")
    parser.add_argument("--cpus", type=int, default=1, help="CPU to use")

    args = parser.parse_args()

    if args.cpus < 1:
        raise ValueError("CPUS must be greater than zero!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .tar!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    tar_files: typing.List[str] = list()

    input_data = step00.read_pickle(args.input).T
    taxa = list(input_data.columns)
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    input_data["Premature"] = list(map(lambda x: metadata.loc[x, "Premature"], input_data.index))
    orders = ("PTB", "Normal")
    print(input_data)

    classifier = sklearn.ensemble.RandomForestClassifier(max_features=None, n_jobs=args.cpus, random_state=0, verbose=1)

    for site in set(input_data["Site"]):
        tmp_data = input_data.loc[(input_data["Site"] == site)]

        # Get Feature Importances
        classifier.fit(tmp_data[taxa], tmp_data["Premature"])
        feature_importances = classifier.feature_importances_
        best_features = list(map(lambda x: x[1], sorted(zip(feature_importances, taxa), reverse=True)))

        # Draw Feature Importances
        matplotlib.use("Agg")
        matplotlib.rcParams.update(step00.matplotlib_parameters)
        seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
        seaborn.distplot(list(filter(lambda x: x > 0, feature_importances)), hist=True, kde=False, rug=True, ax=ax)

        matplotlib.pyplot.title("Feature Importances on {0}".format(site))
        matplotlib.pyplot.xlabel("Feature Importances")
        matplotlib.pyplot.ylabel("Number of Features")
        matplotlib.pyplot.grid(True)
        tar_files.append("{0}+importances.pdf".format(site))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

        # Select best features
        flag = True
        while flag:
            print(len(best_features), "features!!")

            flag = False
            classifier.fit(tmp_data[best_features], tmp_data["Premature"])
            feature_importances = list(classifier.feature_importances_)
            best_features = list(map(lambda x: x[1], sorted(filter(lambda x: x[0] > 0, zip(feature_importances, best_features)), reverse=True)))

            if list(filter(lambda x: x == 0, feature_importances)):
                flag = True
            else:
                tar_files.append("{0}+importances.tsv".format(site))
                with open(tar_files[-1], "w") as f:
                    f.write("Feature\tImportance\n")
                    for feature, importance in sorted(filter(lambda x: x[0] > 0, zip(feature_importances, best_features)), reverse=True):
                        f.write("{0}\t{1}\t".format(feature, importance))

        # Run K-fold
        k_fold = sklearn.model_selection.StratifiedKFold(n_splits=5)
        test_scores = list()
        for i in range(1, 101):
            print("With", i, "/ 100% !!")
            used_columns = taxa[:len(taxa) * i // 100]
            for j, (train_index, test_index) in enumerate(k_fold.split(tmp_data[used_columns], tmp_data["Premature"])):
                x_train, x_test = tmp_data.iloc[train_index][used_columns], tmp_data.iloc[test_index][used_columns]
                y_train, y_test = tmp_data.iloc[train_index]["Premature"], tmp_data.iloc[test_index]["Premature"]

                classifier.fit(x_train, y_train)

                for metric in step00.derivations:
                    test_scores.append((i, metric, step00.aggregate_confusion_matrix(sklearn.metrics.confusion_matrix(y_test, classifier.predict(x_test)), metric)))

        # Draw K-fold
        score_data = pandas.DataFrame.from_records(test_scores, columns=["Features", "Metrics", "Values"])
        matplotlib.use("Agg")
        matplotlib.rcParams.update(step00.matplotlib_parameters)
        seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
        seaborn.lineplot(data=score_data, x="Features", y="Values", hue="Metrics", style="Metrics", ax=ax, markers=True, markersize=20)
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.ylim(0, 1)
        ax.invert_xaxis()
        tar_files.append("{0}+metrics.pdf".format(site))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

        for i, feature in enumerate(taxa[:10]):
            print(i, feature)
            matplotlib.use("Agg")
            matplotlib.rcParams.update(step00.matplotlib_parameters)
            seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
            seaborn.violinplot(data=tmp_data, x="Premature", y=feature, order=orders, ax=ax, inner="box")
            statannot.add_stat_annotation(ax, data=tmp_data, x="Premature", y=feature, order=orders, box_pairs=itertools.combinations(orders, 2), text_format="simple", loc="inside", verbose=2, test="t-test_ind")

            matplotlib.pyplot.ylabel(step00.simplified_taxonomy(feature))

            tar_files.append("{0}+Violin_{1}.pdf".format(site, i))
            fig.savefig(tar_files[-1])
            matplotlib.pyplot.close(fig)

    # Save data
    with tarfile.open(args.output, "w") as tar:
        for file_name in tar_files:
            print(file_name)
            tar.add(file_name, arcname=file_name)
