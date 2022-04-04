"""
step26-1.py: RandomForest Classifier
"""
import argparse
import itertools
import tarfile
import typing
import matplotlib
import matplotlib.pyplot
import numpy
import seaborn
import sklearn.ensemble
import sklearn.manifold
import sklearn.model_selection
import sklearn.tree
import statannotations.Annotator
import pandas
import tqdm
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Train TAR.gz file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output basename")
    parser.add_argument("--cpus", type=int, default=1, help="CPUs to use")
    parser.add_argument("--split", type=int, default=5, help="K-fold split")

    args = parser.parse_args()

    if args.cpus < 1:
        raise ValueError("CPUS must be greater than zero!!")
    elif not args.metadata.endswith(".tsv"):
        raise ValueError("Metadata file must end with .TSV!!")
    elif not args.output.endswith(".tar"):
        raise ValueError("Output file must end with .tar!!")
    elif args.split < 2:
        raise ValueError("SPLIT must be greater than one!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    tar_files: typing.List[str] = list()

    input_data = step00.read_pickle(args.input).T
    taxa = list(input_data.columns)
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    target = "Simple Premature"
    orders = ["Early PTB", "Late PTB+Normal"]
    # target = "Premature"
    # orders = ["PTB", "Normal"]

    input_data = pandas.concat([input_data, metadata], axis="columns", join="inner", verify_integrity=True)
    print(input_data)

    classifier = sklearn.ensemble.RandomForestClassifier(max_features=None, n_jobs=args.cpus, random_state=0, verbose=1)
    k_fold = sklearn.model_selection.StratifiedKFold(n_splits=args.split)

    for site in tqdm.tqdm(step00.selected_long_sites):
        tmp_data = input_data.loc[(input_data["Site"] == site)]

        if len(tmp_data) < args.split:
            continue
        elif len(set(tmp_data[target])) != len(orders):
            continue

        # Get Feature Importances
        classifier.fit(tmp_data[taxa], tmp_data[target])
        feature_importances = list(classifier.feature_importances_)
        best_features = list(map(lambda x: x[1], sorted(zip(feature_importances, taxa), reverse=True)))

        # Draw Feature Importances
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
        seaborn.distplot(feature_importances, hist=True, kde=False, rug=True, ax=ax)

        matplotlib.pyplot.title("Feature Importances on {0}".format(site))
        matplotlib.pyplot.xlabel("Feature Importances")
        matplotlib.pyplot.ylabel("Number of Features")
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.tight_layout()

        tar_files.append("{0}+importances.pdf".format(site))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

        # Select best features
        test_scores = list()
        flag = True

        for j, (train_index, test_index) in enumerate(k_fold.split(tmp_data[best_features], tmp_data[target])):
            x_train, x_test = tmp_data.iloc[train_index][best_features], tmp_data.iloc[test_index][best_features]
            y_train, y_test = tmp_data.iloc[train_index][target], tmp_data.iloc[test_index][target]

            classifier.fit(x_train, y_train)

            for metric in step00.derivations:
                try:
                    test_scores.append((len(best_features), metric, step00.aggregate_confusion_matrix(numpy.sum(sklearn.metrics.multilabel_confusion_matrix(y_test, classifier.predict(x_test)), axis=0), metric)))
                except AssertionError:
                    pass

        while flag:
            print(len(best_features), "features!!")

            flag = False
            classifier.fit(tmp_data[best_features], tmp_data[target])
            feature_importances = list(classifier.feature_importances_)
            best_features = list(map(lambda x: x[1], sorted(filter(lambda x: x[0] > 0, zip(feature_importances, best_features)), reverse=True)))

            if not best_features:
                break

            for j, (train_index, test_index) in enumerate(k_fold.split(tmp_data[best_features], tmp_data[target])):
                x_train, x_test = tmp_data.iloc[train_index][best_features], tmp_data.iloc[test_index][best_features]
                y_train, y_test = tmp_data.iloc[train_index][target], tmp_data.iloc[test_index][target]

                classifier.fit(x_train, y_train)

                for metric in step00.derivations:
                    try:
                        test_scores.append((len(best_features), metric, step00.aggregate_confusion_matrix(sklearn.metrics.confusion_matrix(y_test, classifier.predict(x_test)), metric)))
                    except AssertionError:
                        continue

            if list(filter(lambda x: x == 0, feature_importances)):
                flag = True

        best_BA, tmp_features = -1.0, best_features[:]

        for i in range(1, len(best_features)):
            for j, (train_index, test_index) in enumerate(k_fold.split(tmp_data[best_features[:i]], tmp_data[target])):
                x_train, x_test = tmp_data.iloc[train_index][best_features[:i]], tmp_data.iloc[test_index][best_features[:i]]
                y_train, y_test = tmp_data.iloc[train_index][target], tmp_data.iloc[test_index][target]

                classifier.fit(x_train, y_train)

                for metric in step00.derivations:
                    try:
                        test_scores.append((i, metric, step00.aggregate_confusion_matrix(numpy.sum(sklearn.metrics.multilabel_confusion_matrix(y_test, classifier.predict(x_test)), axis=0), metric)))
                        if metric == "Balanced Accuracy":
                            BA = step00.aggregate_confusion_matrix(numpy.sum(sklearn.metrics.multilabel_confusion_matrix(y_test, classifier.predict(x_test)), axis=0), metric)
                            if best_BA < BA:
                                best_BA = BA
                                tmp_features = best_features[:]
                    except AssertionError:
                        continue

        heatmap_data = pandas.DataFrame(data=numpy.zeros((len(orders), len(orders))), index=orders, columns=orders, dtype=int)

        for j, (train_index, test_index) in enumerate(k_fold.split(tmp_data[tmp_features], tmp_data[target])):
            x_train, x_test = tmp_data.iloc[train_index][tmp_features], tmp_data.iloc[test_index][tmp_features]
            y_train, y_test = tmp_data.iloc[train_index][target], tmp_data.iloc[test_index][target]

            classifier.fit(x_train, y_train)

            for real, prediction in zip(y_test, classifier.predict(x_test)):
                heatmap_data.loc[real, prediction] += 1

        fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))

        seaborn.heatmap(data=heatmap_data, annot=True, fmt="d", cbar=False, square=True, xticklabels=True, yticklabels=True, ax=ax)

        # Draw K-fold
        score_data = pandas.DataFrame.from_records(test_scores, columns=["Features", "Metrics", "Values"])
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

        seaborn.lineplot(data=score_data, x="Features", y="Values", hue="Metrics", style="Metrics", ax=ax, markers=True, markersize=20)
        matplotlib.pyplot.axvline(x=len(tmp_features), linestyle="--", color="k")
        matplotlib.pyplot.text(x=len(tmp_features), y=0.2, s=f"Best BA {best_BA} with {len(tmp_features)} features", fontsize="xx-small", color="k", horizontalalignment="center", verticalalignment="bottom", rotation="vertical")

        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.ylim(0, 1)
        matplotlib.pyplot.title(site)
        ax.invert_xaxis()
        matplotlib.pyplot.tight_layout()

        tar_files.append("{0}+metrics.pdf".format(site))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

        for i, feature in enumerate(taxa[:10]):
            fig, ax = matplotlib.pyplot.subplots(figsize=(24, 24))
            seaborn.violinplot(data=tmp_data, x=target, y=feature, order=orders, ax=ax, inner="box")
            try:
                statannotations.Annotator.Annotator(ax, list(itertools.combinations(orders, 2)), data=tmp_data, x=target, y=feature, order=orders).configure(test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()
            except ValueError:
                pass

            matplotlib.pyplot.ylabel(step00.consistency_taxonomy(feature, 1))
            matplotlib.pyplot.title(site)
            matplotlib.pyplot.tight_layout()

            tar_files.append("{0}+Violin_{1}.pdf".format(site, i))
            fig.savefig(tar_files[-1])
            matplotlib.pyplot.close(fig)

    # Save data
    with tarfile.open(args.output, "w") as tar:
        for file_name in tqdm.tqdm(tar_files):
            tar.add(file_name, arcname=file_name)
