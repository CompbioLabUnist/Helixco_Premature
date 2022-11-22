"""
step62-1.py: RandomForest Classifier with DAT
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

ratio_threshold = 2
p_threshold = 0.05

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=str, help="Train TSV file")
    parser.add_argument("DAT", type=str, help="DAT TSV file")
    parser.add_argument("metadata", type=str, help="Metadata TSV file")
    parser.add_argument("output", type=str, help="Output basename")
    parser.add_argument("--cpus", type=int, default=1, help="CPUs to use")
    parser.add_argument("--split", type=int, default=5, help="K-fold split")

    args = parser.parse_args()

    if args.cpus < 1:
        raise ValueError("CPUS must be greater than zero!!")
    elif not args.input.endswith(".tsv"):
        raise ValueError("INPUT must end with .TSV!!")
    elif not args.DAT.endswith(".tsv"):
        raise ValueError("DAT must end with .TSV!!")
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

    DAT_data = pandas.read_csv(args.DAT, sep="\t", index_col=0)
    DAT_data = DAT_data.loc[((DAT_data["log2FoldChange"] > numpy.log2(ratio_threshold)) | (DAT_data["log2FoldChange"] < -1 * numpy.log2(ratio_threshold))) & (DAT_data["padj"] < p_threshold)]
    print(DAT_data)

    input_data = pandas.read_csv(args.input, sep="\t", index_col=0).T
    input_data = input_data.loc[:, DAT_data.index]
    input_data = input_data.loc[:, list(filter(step00.filtering_taxonomy, list(input_data.columns)))]
    taxa = list(input_data.columns)
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    input_data = pandas.concat([input_data, metadata], axis="columns", join="inner", verify_integrity=True)
    print(input_data)

    target = "Simple Premature"
    orders = ["Early PTB", "Late PTB+Normal"]

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
    for i in tqdm.trange(1, len(best_features)):
        for j, (train_index, test_index) in enumerate(k_fold.split(tmp_data[best_features[:i]], tmp_data[target])):
            x_train, x_test = tmp_data.iloc[train_index][best_features[:i]], tmp_data.iloc[test_index][best_features[:i]]
            y_train, y_test = tmp_data.iloc[train_index][target], tmp_data.iloc[test_index][target]

            classifier.fit(x_train, y_train)

            for metric in step00.derivations:
                try:
                    test_scores.append((i, metric, step00.aggregate_confusion_matrix(numpy.sum(sklearn.metrics.multilabel_confusion_matrix(y_test, classifier.predict(x_test)), axis=0), metric)))
                except ZeroDivisionError:
                    continue

    score_data = pandas.DataFrame.from_records(test_scores, columns=["Features", "Metrics", "Values"])
    best_BA, tmp_features = -1.0, best_features[:]

    for i in sorted(set(score_data["Features"])):
        BA = numpy.mean(score_data.loc[(score_data["Features"] == i) & (score_data["Metrics"] == "Balanced Accuracy"), "Values"])
        if best_BA < BA:
            best_BA = BA
            tmp_features = best_features[:i]

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

    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.ylabel("Evaluations")
    ax.invert_xaxis()
    matplotlib.pyplot.tight_layout()

    tar_files.append("metrics.pdf")
    fig.savefig(tar_files[-1])
    matplotlib.pyplot.close(fig)

    for derivation in step00.derivations:
        print("--", derivation, numpy.mean(score_data.loc[(score_data["Features"] == len(tmp_features)) & (score_data["Metrics"] == derivation), "Values"]), numpy.std(score_data.loc[(score_data["Features"] == len(tmp_features)) & (score_data["Metrics"] == derivation), "Values"]))

    for i, feature in enumerate(best_features[:10]):
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
