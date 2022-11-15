"""
step56-1.py: K-Nearest Neighbors Classifier
"""
import argparse
import itertools
import multiprocessing
import tarfile
import typing
import imblearn.over_sampling
import matplotlib
import matplotlib.pyplot
import numpy
import seaborn
import sklearn.manifold
import sklearn.model_selection
import sklearn.neighbors
import sklearn.tree
import statannotations.Annotator
import pandas
import tqdm
import step00

tmp_data = pandas.DataFrame()


def std(column: str) -> float:
    return numpy.std(tmp_data[column])


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

    input_data = pandas.read_csv(args.input, sep="\t", skiprows=1, index_col=["#OTU ID"]).groupby("taxonomy").sum().T
    input_data = input_data.loc[:, list(filter(step00.filtering_taxonomy, list(input_data.columns)))]
    taxa = list(input_data.columns)
    print(input_data)

    metadata = pandas.read_csv(args.metadata, sep="\t", skiprows=[1]).dropna(axis="columns", how="all").set_index(keys="#SampleID", verify_integrity=True)
    print(metadata)

    input_data = pandas.concat([input_data, metadata], axis="columns", join="inner", verify_integrity=True)
    print(input_data)

    target = "Simple Premature"
    orders = ["Early PTB", "Late PTB+Normal"]
    # target = "Premature"
    # orders = ["PTB", "Normal"]

    classifier = sklearn.neighbors.KNeighborsClassifier(algorithm="brute", weights="distance", n_jobs=args.cpus)
    k_fold = sklearn.model_selection.StratifiedKFold(n_splits=args.split)
    oversampler = imblearn.over_sampling.RandomOverSampler(random_state=0)

    for site in tqdm.tqdm(step00.selected_long_sites):
        selected_data = input_data.loc[(input_data["Site"] == site)]

        X, y = oversampler.fit_resample(selected_data[taxa].to_numpy(), selected_data[target])
        tmp_data = pandas.DataFrame(X, columns=taxa)
        tmp_data[target] = y

        if len(tmp_data) < args.split:
            continue
        elif len(set(tmp_data[target])) != len(orders):
            continue

        flag = False
        for order in orders:
            if len(tmp_data.loc[(tmp_data[target] == order)]) < args.split:
                flag = True
        if flag:
            continue

        # Get Feature Importances
        with multiprocessing.Pool(args.cpus) as pool:
            feature_stds = pool.map(std, taxa)
        best_features = list(map(lambda x: x[1], sorted(zip(feature_stds, taxa), reverse=True)))

        # Draw Feature STDs
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))
        seaborn.distplot(feature_stds, hist=True, kde=False, rug=True, ax=ax)

        matplotlib.pyplot.title("Feature STDs on {0}".format(site))
        matplotlib.pyplot.xlabel("Feature Standard Deviations")
        matplotlib.pyplot.ylabel("Number of Features")
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.tight_layout()

        tar_files.append("{0}+STDs.pdf".format(site))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

        # Select best features
        test_scores = list()
        for i in tqdm.trange(1, 200):
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
        matplotlib.pyplot.title(site)
        matplotlib.pyplot.tight_layout()
        tar_files.append("{0}+heatmap.pdf".format(site))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

        # Draw K-fold
        fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

        seaborn.lineplot(data=score_data, x="Features", y="Values", hue="Metrics", style="Metrics", ax=ax, markers=True, markersize=20)
        matplotlib.pyplot.axvline(x=len(tmp_features), linestyle="--", color="k")
        matplotlib.pyplot.text(x=len(tmp_features), y=0.1, s=f"Best BA {best_BA:.3f} with {len(tmp_features)} features", fontsize="xx-small", color="k", horizontalalignment="right", verticalalignment="baseline", rotation="vertical")

        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.ylim(0, 1)
        matplotlib.pyplot.title(site)
        matplotlib.pyplot.ylabel("Evaluations")
        ax.invert_xaxis()
        matplotlib.pyplot.tight_layout()

        tar_files.append("{0}+metrics.pdf".format(site))
        fig.savefig(tar_files[-1])
        matplotlib.pyplot.close(fig)

        for i, feature in enumerate(best_features[:10]):
            print(feature)

            fig, ax = matplotlib.pyplot.subplots(figsize=(18, 18))
            seaborn.violinplot(data=input_data, x=target, y=feature, order=orders, ax=ax, inner="box", cut=1, linewidth=5)
            try:
                statannotations.Annotator.Annotator(ax, list(itertools.combinations(orders, 2)), data=input_data, x=target, y=feature, order=orders).configure(test="Mann-Whitney", text_format="simple", loc="inside", verbose=0).apply_and_annotate()
            except ValueError:
                pass

            matplotlib.pyplot.ylabel(f"{step00.simplified_taxonomy(feature)} in {site}")
            matplotlib.pyplot.tight_layout()

            tar_files.append("{0}+Violin_{1}.pdf".format(site, i))
            fig.savefig(tar_files[-1])
            matplotlib.pyplot.close(fig)

    # Save data
    with tarfile.open(args.output, "w") as tar:
        for file_name in tqdm.tqdm(tar_files):
            tar.add(file_name, arcname=file_name)
