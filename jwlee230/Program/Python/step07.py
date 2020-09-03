"""
step07.py: make metadata
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("info", help="Informative XLSX file", type=str, nargs=1)
    parser.add_argument("input", help="Input file (FQ.gz) names", type=str, nargs="+")

    args = parser.parse_args()

    if not args.info[0].endswith(".xlsx"):
        raise ValueError("INFO file must end with .xlsx")

    info_data = pandas.read_excel(args.info[0])

    data = pandas.DataFrame()
    data["#SampleID"] = sorted(list(set(list(map(lambda x: x.split("_")[0], list(map(lambda x: x.split("/")[-1], args.input)))))))
    data["BarcodeSequence"] = ""
    data["LinkPrimerSequence"] = ""
    data["Site"] = list(map(lambda x: {"B1": "Baby-1day", "B3": "Baby-3day", "B5": "Baby-5day", "M": "Mouth", "C": "Cervix", "V": "Vagina", "P": "Placenta"}[x.split("-")[-1]], data["#SampleID"]))
    data["Mother"] = list(map(lambda x: x.split("-")[0], data["#SampleID"]))
    data["BabyNumber"] = list(map(lambda x: x.split("-")[1] if len(x.split("-")) == 3 else "1", data["#SampleID"]))
    data["Premature"] = list(map(lambda x: "Premature" if int(x) < 37 else "Normal", list(map(lambda x: list(info_data.loc[(info_data["바코드"].str.startswith(x))]["분만주수"])[0].split("+")[0], data["Mother"]))))
    data["Description"] = ""

    print("\t".join(data.columns))
    print("#q2:types", "\t".join(["categorical"] * (len(data.columns) - 1)), sep="\t")
    for index, row in data.iterrows():
        print("\t".join(row))
