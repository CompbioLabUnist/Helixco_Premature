"""
step67-2.py: RandomForest Classifier with DAT-union
"""
import argparse
import itertools
import tarfile
import typing
import imblearn.over_sampling
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

    parser.add_argument("input", type=str, help="Train TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output basename")
    parser.add_argument("--cpus", type=int, default=1, help="CPUs to use")
    parser.add_argument("--split", type=int, default=5, help="K-fold split")

    args = parser.parse_args()

    if args.cpus < 1:
        raise ValueError("CPUS must be greater than zero!!")
    elif not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
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

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    input_data = input_data.loc[:, list(filter(step00.filtering_taxonomy, list(input_data.columns)))]
    print(input_data)

    if input_data.empty:
        exit()

    taxa = list(input_data.columns)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    input_data = pandas.concat([input_data, metadata], axis="columns", join="inner", verify_integrity=True)
    print(input_data)

    target = "Premature"
    orders = ["PTB", "Normal"]

    classifier = sklearn.ensemble.RandomForestClassifier(max_features=None, n_jobs=args.cpus, random_state=0, verbose=1)
    k_fold = sklearn.model_selection.StratifiedKFold(n_splits=args.split)
    oversampler = imblearn.over_sampling.SMOTE(random_state=42, k_neighbors=3, n_jobs=args.cpus)

    X, y = oversampler.fit_resample(input_data[taxa].to_numpy(), input_data[target])
    tmp_data = pandas.DataFrame(X, columns=taxa)
    tmp_data[target] = y

    classifier.fit(tmp_data[taxa], tmp_data[target])
    feature_importances = list(classifier.feature_importances_)
    best_features = list(map(lambda x: x[1], sorted(zip(feature_importances, taxa), reverse=True)))

    test_scores = list()
    for i in tqdm.trange(1, len(best_features) + 1):
        for j, (train_index, test_index) in enumerate(k_fold.split(tmp_data[best_features[:i]], tmp_data[target])):
            x_train, x_test = tmp_data.iloc[train_index][best_features[:i]], tmp_data.iloc[test_index][best_features[:i]]
            y_train, y_test = tmp_data.iloc[train_index][target], tmp_data.iloc[test_index][target]

            classifier.fit(x_train, y_train)

            for metric in step00.derivations:
                try:
                    test_scores.append((i, metric, step00.aggregate_confusion_matrix(sklearn.metrics.confusion_matrix(y_test, classifier.predict(x_test)), metric)))
                except ZeroDivisionError:
                    continue

    score_data = pandas.DataFrame.from_records(test_scores, columns=["Features", "Metrics", "Values"])
    best_BA, tmp_features = -1.0, best_features[:]

    for i in sorted(set(score_data["Features"])):
        BA = numpy.mean(score_data.loc[(score_data["Features"] == i) & (score_data["Metrics"] == "BA"), "Values"])
        if best_BA < BA:
            best_BA = BA
            tmp_features = best_features[:i]

    # Importances
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

    seaborn.histplot(data=feature_importances, stat="count", kde=True, ax=ax)

    matplotlib.pyplot.xlabel("Importances")
    matplotlib.pyplot.title("PTB vs. Normal")
    matplotlib.pyplot.tight_layout()

    tar_files.append("importances.pdf")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    importance_data = pandas.DataFrame(index=best_features, data=sorted(feature_importances, reverse=True), columns=["Importance"]).sort_index()
    importance_data.index.name = "Taxonomy"
    print(importance_data)
    importance_data.to_csv(args.output.replace(".tar", ".importance.tsv"), sep="\t")

    # Heatmap
    heatmap_data = pandas.DataFrame(data=numpy.zeros((len(orders), len(orders))), index=orders, columns=orders, dtype=int)
    for j, (train_index, test_index) in enumerate(k_fold.split(tmp_data[tmp_features], tmp_data[target])):
        x_train, x_test = tmp_data.iloc[train_index][tmp_features], tmp_data.iloc[test_index][tmp_features]
        y_train, y_test = tmp_data.iloc[train_index][target], tmp_data.iloc[test_index][target]

        classifier.fit(x_train, y_train)

        for real, prediction in zip(y_test, classifier.predict(x_test)):
            heatmap_data.loc[real, prediction] += 1

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.heatmap(data=heatmap_data, annot=True, fmt="d", cbar=False, square=True, xticklabels=True, yticklabels=True, ax=ax)

    matplotlib.pyplot.xlabel("Prediction")
    matplotlib.pyplot.ylabel("Real")
    matplotlib.pyplot.tight_layout()
    tar_files.append("heatmap.pdf")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    # Draw K-fold
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

    seaborn.lineplot(data=score_data, x="Features", y="Values", hue="Metrics", style="Metrics", ax=ax, markers=True, markersize=20)
    matplotlib.pyplot.axvline(x=len(tmp_features), linestyle="--", color="k")
    matplotlib.pyplot.text(x=len(tmp_features), y=0.1, s=f"Best BA {best_BA:.3f} with {len(tmp_features)} features", fontsize="xx-small", color="k", horizontalalignment="right", verticalalignment="baseline", rotation="vertical")

    matplotlib.pyplot.xticks(range(1, len(taxa) + 1), range(1, len(taxa) + 1))
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.ylabel("Evaluations")
    matplotlib.pyplot.title("PTB vs. Normal")
    ax.invert_xaxis()
    matplotlib.pyplot.tight_layout()

    tar_files.append("metrics.pdf")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    # Draw bar
    tmp_data = score_data.loc[(score_data["Features"] == len(tmp_features))]

    fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))

    seaborn.barplot(data=tmp_data, x="Metrics", y="Values", order=step00.derivations, ax=ax)

    matplotlib.pyplot.xlabel("")
    matplotlib.pyplot.ylabel("Evaluations")
    matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.title("PTB vs. Normal")
    matplotlib.pyplot.tight_layout()

    tar_files.append("bar.pdf")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    raw_evaluation_data: typing.List[typing.List[str]] = list()
    for i in range(1, len(best_features) + 1):
        tmp = list()
        for derivation in step00.derivations:
            d = score_data.loc[(score_data["Features"] == i) & (score_data["Metrics"] == derivation)]
            tmp.append(f"{numpy.mean(d)}±{numpy.std(d)}")
        raw_evaluation_data.append(tmp)
    evaluation_data = pandas.DataFrame(raw_evaluation_data, columns=step00.derivations)
    print(evaluation_data)
    evaluation_data.to_csv(args.output.replace(".tar", ".evaluation.tsv"), sep="\t")

    for derivation in step00.derivations:
        print("--", derivation, numpy.mean(score_data.loc[(score_data["Features"] == len(tmp_features)) & (score_data["Metrics"] == derivation), "Values"]), numpy.std(score_data.loc[(score_data["Features"] == len(tmp_features)) & (score_data["Metrics"] == derivation), "Values"]))

    for i, feature in enumerate(best_features):
        print(feature)

        fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
        seaborn.violinplot(data=input_data, x=target, y=feature, order=orders, ax=ax, inner="box", cut=1, linewidth=5)
        try:
            statannotations.Annotator.Annotator(ax, list(itertools.combinations(orders, 2)), data=input_data, x=target, y=feature, order=orders).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
        except ValueError:
            pass

        matplotlib.pyplot.ylabel(f"{step00.simplified_taxonomy(feature)}")
        matplotlib.pyplot.tight_layout()

        tar_files.append(f"Violin_{i}.pdf")
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

    # Save data
    with tarfile.open(args.output, "w") as tar:
        for file_name in tqdm.tqdm(tar_files):
            tar.add(file_name, arcname=file_name)
