"""
step07.py: make metadata
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input file (FQ.gz) names", type=str, nargs="+")

    args = parser.parse_args()

    data = pandas.DataFrame()
    data["#SampleID"] = sorted(list(set(list(map(lambda x: x.split("_")[0], list(map(lambda x: x.split("/")[-1], args.input)))))))
    data["BarcodeSequence"] = ""
    data["LinkPrimerSequence"] = ""
    data["Site"] = list(map(lambda x: {"B1": "Baby-1day", "B3": "Baby-3day", "B5": "Baby-5day", "M": "Mouth", "C": "Cervix", "V": "Vagina", "P": "Placenta"}[x.split("-")[-1]], data["#SampleID"]))
    data["Mother"] = list(map(lambda x: x.split("-")[0], data["#SampleID"]))
    data["BabyNumber"] = list(map(lambda x: x.split("-")[1] if len(x.split("-")) == 3 else "1", data["#SampleID"]))
    data["Description"] = ""

    print("\t".join(data.columns))
    print("#q2:types", "\t".join(["categorical"] * (len(data.columns) - 1)), sep="\t")
    for index, row in data.iterrows():
        print("\t".join(row))
