"""
step07-5.py: make metadata for no-twin
"""
import argparse
import pandas
import step00

metadata = pandas.DataFrame()


def change_index(ID: str) -> str:
    a, b, c, d = ID.split("-")
    if a == "Stool":
        a = "First"

    if (a == "First") and (c == "Mother"):
        if len(metadata.loc[(metadata["Data"] == "First") & (metadata["Mother"] == b)].to_numpy()) == 1:
            c = "0"
        else:
            c = "1"
    elif (a in ["Second", "Third"]) and (c == "Mother"):
        if len(metadata.loc[(metadata["Data"].isin(["Second", "Third"])) & (metadata["Mother"] == b)].to_numpy()) == 1:
            c = "0"
        else:
            c = "1"

    if a == "First":
        return "-".join((a, b, c))
    elif a in ["Second", "Third"]:
        return "-".join(metadata.loc[(metadata["Mother"] == b) & (metadata["Neonate"] == c), ["Data", "Mother", "Neonate"]].to_numpy()[0])
    else:
        raise Exception("Something went wrong!!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input TSV file", type=str)
    parser.add_argument("info", help="Informative XLSX file", type=str)

    args = parser.parse_args()

    if not args.input.endswith(".tsv"):
        raise ValueError("Input file must end with .TSV!!")
    elif not args.info.endswith(".xlsx"):
        raise ValueError("INFO file must end with .XLSX!!")

    input_data = pandas.read_csv(args.input, sep="\t", comment="#")

    metadata = pandas.read_excel(args.info, sheet_name=0, dtype=str)
    metadata["ID"] = list(map(lambda x: "-".join(x), zip(metadata["Data"], metadata["Mother"], metadata["Neonate"])))
    metadata.set_index("ID", inplace=True, verify_integrity=True)

    singleton_data = metadata.loc[(metadata["Neonate"] == "0")]
    singleton_mother = set(map(lambda x: "-".join(x), zip(singleton_data["Data"], singleton_data["Mother"])))

    data = pandas.DataFrame()
    data["#SampleID"] = sorted(filter(lambda x: "-".join(x.split("-")[:2]) in singleton_mother, input_data["sample-id"]))
    data["BarcodeSequence"] = ""
    data["LinkPrimerSequence"] = ""
    data["Site"] = list(map(lambda x: {"B1": "Neonate-1day", "B3": "Neonate-3day", "B5": "Neonate-5day", "M": "Mouth", "C": "Cervix", "V": "Vagina", "P": "Placenta", "S1": "Stool-1day", "S3": "Stool-3day", "S5": "Stool-5day"}[x.split("-")[-1]], data["#SampleID"]))
    data["ID"] = list(map(change_index, data["#SampleID"]))
    for c in list(metadata.columns):
        data[c] = list(map(lambda x: metadata.loc[x, c], data["ID"]))
    data["Simple Premature"] = list(map(lambda x: "Early PTB" if (x == "Early PTB") else "Late PTB+Normal", data["Detail Premature"]))
    data["Description"] = ""
    del data["ID"]

    print("\t".join(list(data.columns)))
    print("#q2:types", "\t".join(list(map(lambda x: "numeric" if x in step00.numeric_columns else "categorical", list(data.columns)[1:]))), sep="\t")
    for index, row in data.iterrows():
        print("\t".join(list(map(str, row))))
