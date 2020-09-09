"""
step25.py: make pseudo-sample
"""
import argparse
import random
import pandas
import step00

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TAR.gz file", type=str, nargs=1)
    parser.add_argument("output", help="Output TAR.gz file", type=str, nargs=1)
    parser.add_argument("--dev", help="deviation for normal distribution", type=float, default=0.05)

    args = parser.parse_args()

    data: pandas.DataFrame = step00.read_pickle(args.input[0])

    random.seed(a=0)
    data = data.applymap(lambda x: x * (1 + random.normalvariate(0, args.dev) if x > 0 else abs(random.normalvariate(0, args.dev))))

    step00.make_pickle(args.output[0], data)
