"""
step38-1.py: Violin plots with sites for Macrogen data
"""
import argparse
import itertools
import multiprocessing
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import scipy.stats
import seaborn
import step00


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input tar.gz file", type=str)
    parser.add_argument("output", help="Output PDF file", type=str)
    parser.add_argument("--cpus", help="Number of cpus", type=int, default=1)

    args = parser.parse_args()

    if not args.input.endswith(".tar.gz"):
        raise ValueError("Input file must end with .tar.gz!!")
    elif not args.output.endswith(".pdf"):
        raise ValueError("Output file must end with .PDF!!")
    elif args.cpus < 1:
        raise ValueError("CPUS must be a positive integer!!")

    matplotlib.use("Agg")
    matplotlib.rcParams.update(step00.matplotlib_parameters)
    seaborn.set_theme(context="poster", style="whitegrid", rc=step00.matplotlib_parameters)

    input_data = step00.read_pickle(args.input)
    input_data.index = list(map(step00.simplified_taxonomy, list(input_data.index)))
    input_data = input_data.iloc[:, 1:].T
    print(input_data)

    g = seaborn.clustermap(data=input_data, z_score=1, figsize=(32, 18), row_cluster=True, col_cluster=True, xticklabels=False, yticklabels=False, cmap="bwr", vmin=-3, center=0, vmax=3)
    g.ax_heatmap.set_xlabel("Taxonomy")
    g.ax_heatmap.set_ylabel("Sample")

    g.savefig(args.output)
