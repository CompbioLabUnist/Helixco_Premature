"""
step07-2.py: make metadata
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("info", help="Informative TSV file", type=str)
    parser.add_argument("input", help="Input FASTQ.gz files", type=str, nargs="+")

    args = parser.parse_args()

    if not args.info.endswith(".tsv"):
        raise ValueError("INFO file must end with .tsv")

    info_data = pandas.read_csv(args.info, sep="\t")

    data = pandas.DataFrame()
    data["#SampleID"] = info_data["run_accession"]
    data["BarcodeSequence"] = ""
    data["LinkPrimerSequence"] = ""
    data["site"] = info_data["scientific_name"]
    data["Description"] = ""

    print("\t".join(data.columns))
    print("#q2:types", "\t".join(["categorical"] * (len(data.columns) - 1)), sep="\t")
    for index, row in data.iterrows():
        print("\t".join(row))
