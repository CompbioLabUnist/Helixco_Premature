"""
step01-4.py: make manifest file from stool data
"""
import sys

f = sys.argv[1:]
f1 = sorted(list(filter(lambda x: x.endswith("R1_001.fastq.gz"), f)))
f2 = sorted(list(filter(lambda x: x.endswith("R2_001.fastq.gz"), f)))
assert len(f1) == len(f2)

print("sample-id\tforward-absolute-filepath\treverse-absolute-filepath")
for fa, fb in zip(f1, f2):
    name = fa[:-19]
    if fb.startswith(name):
        print(name[name.rfind("/") + 1:].replace("_", "-"), fa, fb, sep="\t")
    else:
        raise ValueError(fa, fb)
