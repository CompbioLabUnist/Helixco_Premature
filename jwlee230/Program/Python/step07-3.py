"""
step07.py: make metadata
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("info", help="Informative XLSX file", type=str)
    parser.add_argument("input", help="Input file (Fastq.gz) names", type=str, nargs="+")

    args = parser.parse_args()

    if not args.info.endswith(".xlsx"):
        raise ValueError("INFO file must end with .xlsx")
    elif list(filter(lambda x: not x.endswith(".fastq.gz"), args.input)):
        raise ValueError("INFO file must end with .fastq.gz")

    mother_data = pandas.read_excel(args.info, sheet_name=0)

    numeric_columns = {"Gestational Week", "Weight", "Mother Age", "Hospitalized Day", "Apgar Score", "Weight gain"}

    data = pandas.DataFrame()
    data["#SampleID"] = sorted(list(set(list(map(lambda x: x.split("_")[0], list(map(lambda x: x.split("/")[-1], args.input)))))))
    data["BarcodeSequence"] = ""
    data["LinkPrimerSequence"] = ""
    data["Site"] = list(map(lambda x: {"B1": "Baby-1day", "B3": "Baby-3day", "B5": "Baby-5day", "M": "Mouth", "C": "Cervix", "V": "Vagina", "P": "Placenta"}[x.split("-")[0]], data["#SampleID"]))
    data["Mother"] = list(map(lambda x: "%03d" % int(x.split("-")[1]), data["#SampleID"]))
    data["Baby_Number"] = list(map(lambda x: x.split("-")[1] if len(x.split("-")) == 3 else "1", data["#SampleID"]))

    print(list(data["Mother"]))
    print(list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["분만주수"]), data["Mother"])))

    data["Premature"] = list(map(lambda x: "Premature" if int(x) < 37 else "Normal", list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["분만주수"])[0].split("+")[0], data["Mother"]))))
    data["Detail_Premature"] = list(map(lambda x: "NonlatePremature" if int(x) < 34 else ("LatePremature" if int(x) < 37 else "Normal"), list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["분만주수"])[0].split("+")[0], data["Mother"]))))
    data["Utero_Week"] = list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["분만주수"])[0].split("+")[0], data["Mother"]))
    data["Mother_Age"] = list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["나이"])[0], data["Mother"]))
    data["Mother_Weight"] = list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["체중/키"])[0].split("/")[0], data["Mother"]))
    data["Mother_Height"] = list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["체중/키"])[0].split("/")[1], data["Mother"]))
    data["Mother_BMI"] = list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["BMI"])[0], data["Mother"]))
    data["Early_Labor_Pain"] = list(map(lambda x: "Yes" if "유" in x else "No", list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["조기진통유무"])[0], data["Mother"]))))
    data["PROM"] = list(map(lambda x: "Yes" if "유" in x else "No", list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["조기양막파수 유무"])[0], data["Mother"]))))
    data["C-section"] = list(map(lambda x: "Yes" if "제왕절개" in x else "No", list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["분만형태"])[0], data["Mother"]))))
    data["Steroid"] = list(map(lambda x: "Yes" if "유" in x else "No", list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["산전스테로이드사용여부"])[0], data["Mother"]))))
    data["Antibiotic"] = list(map(lambda x: "Yes" if "유 " in x else "No", list(map(lambda x: list(mother_data.loc[(mother_data["바코드"].str.startswith(x))]["산전항생제(6주이내) 사용여부"])[0], data["Mother"]))))
    data["Description"] = ""

    print("\t".join(data.columns))
    print("#q2:types", "\t".join(list(map(lambda x: "numeric" if x in numeric_columns else "categorical", list(data.columns)[1:]))), sep="\t")
    for index, row in data.iterrows():
        print("\t".join(list(map(str, row))))
