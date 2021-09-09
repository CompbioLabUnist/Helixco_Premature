"""
step09.py: get TSV data from Macrogen excel file
"""
import argparse
import pandas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", help="Input XLSX file(s)", type=str, nargs="+")
    parser.add_argument("metadata", help="Metadata XLSX file", type=str)
    parser.add_argument("output", help="Output TSV file", type=str)

    args = parser.parse_args()

    if list(filter(lambda x: not x.endswith(".xlsx"), args.input)):
        raise ValueError("Input must end with XLSX!!")
    elif not args.metadata.endswith(".xlsx"):
        raise ValueError("Metadata file must end with .xlsx!!")
    elif not args.output.endswith(".tsv"):
        raise ValueError("Output must end with .TSV!!")

    metadata = pandas.read_excel(args.metadata, dtype=str)
    print(metadata)

    data_list = list()
    for input_file in args.input:
        data = pandas.read_excel(input_file, index_col="Organism")

        samples = list(data.columns)[13:]
        data = data.loc[:, samples]
        revise_samples = list()
        print(data)

        for sample in samples:
            sample = sample.split(".")
            name = ""

            if sample[0].isdigit():
                name = "First-"

            if name == "First-":
                if len(sample) == 3:
                    name += "-".join(sample)
                elif (len(sample) == 2) and (sample[1] in ["B1", "B3", "B5"]):
                    name += sample[0] + "-0-" + sample[1]
                elif len(sample) == 2:
                    name += sample[0] + "-Mother-" + sample[1]
                else:
                    raise Exception("Something went wrong!!")
            else:
                tmp_data = metadata.loc[(metadata["Sample"] == sample[0]) & (metadata["ID"] == str(int(sample[1]))), :].to_numpy()[0]
                if int(tmp_data[0]) <= 33:
                    name = "Second-"
                elif (int(tmp_data[0]) == 34) and (sample[0] != "B5"):
                    name = "Second-"
                elif (int(tmp_data[0]) == 35) and (sample[0] != "B3"):
                    name = "Second-"
                else:
                    name = "Third-"

                if tmp_data[1] == "Mother":
                    name += tmp_data[0] + "-Mother-" + tmp_data[2]
                else:
                    name += "-".join(tmp_data[:3])

            revise_samples.append(name)

        data.columns = revise_samples
        data_list.append(data.groupby(by="Organism").sum())

    input_data = pandas.concat(objs=data_list, axis="columns").fillna(value=0)
    input_data["taxonomy"] = list(input_data.index)
    input_data["#Hash"] = list(map(hash, list(input_data.index)))
    input_data.set_index(keys="#Hash", inplace=True, verify_integrity=True)
    print(input_data)

    with open(args.output, "w") as f:
        f.write("# Constructed from xlsx file\n")
        input_data.to_csv(f, sep="\t")
