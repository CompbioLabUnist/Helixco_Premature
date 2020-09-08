"""
step24.py: Draw histogram with parsing
"""
import argparse
import pandas
import matplotlib
import matplotlib.pyplot
import seaborn

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input XLSX file", type=str, nargs=1)
    parser.add_argument("output", help="Output PNG file", type=str, nargs=1)
    parser.add_argument("col", help="Used column", type=str, nargs=1)
    parser.add_argument("num", type=int, help="Which one will be used", choices=[0, 1])
    parser.add_argument("chr", help="How parse with", choices=["/", "+"])

    args = parser.parse_args()

    input_data = pandas.read_excel(args.input[0])

    seaborn.set(context="poster", style="whitegrid")
    fig, ax = matplotlib.pyplot.subplots(figsize=(32, 18))

    seaborn.distplot(list(map(lambda x: x.split(args.chr)[args.num], input_data[args.col[0]])), kde=False, rug=True, axlabel=False)

    fig.savefig(args.output[0])
    matplotlib.pyplot.close(fig)
