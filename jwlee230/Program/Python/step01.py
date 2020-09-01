import sys

f = sys.argv[1:]
f1 = sorted(list(filter(lambda x: x.endswith("1.fastq.gz"), f)))
f2 = sorted(list(filter(lambda x: x.endswith("2.fastq.gz"), f)))

print("sample-id\tforward-absolute-filepath\treverse-absolute-filepath")
for fa, fb in zip(f1, f2):
    name = fa[:-11]
    if fb.startswith(name):
        print(name[name.rfind("/") + 1:], fa, fb, sep="\t")
    else:
        raise ValueError(fa, fb)
