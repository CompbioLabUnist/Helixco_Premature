#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 3)); then
    qiime tools export --input-path $1 --output-path $3
    qiime tools export --input-path $2 --output-path $3
    sed --in-place "1c#OTU ID\ttaxonomy\tconfidence" $3/taxonomy.tsv
    biom add-metadata --input-fp $3/feature-table.biom --observation-metadata-fp $3/taxonomy.tsv --output-fp $3/table.biom --sc-separated "taxonomy"
else
    echo "This script will work if and only if it has **three** arguments"
fi
