"""
step01-2.py: make manifest file
"""
import sys

files = sorted(sys.argv[1:])

print("sample-id\tabsolute-filepath")
for f in files:
    print(f.split("/")[-2], f, sep="\t")
