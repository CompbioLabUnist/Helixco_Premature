#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 3)); then
    qiime tools export --input-path $1 --output-path /tmpfs
    qiime tools export --input-path $2 --output-path /tmpfs
    sed --in-place "1c#OTU ID\ttaxonomy\tconfidence" /tmpfs/taxonomy.tsv
    biom add-metadata --input-fp /tmpfs/feature-table.biom --observation-metadata-fp /tmpfs/taxonomy.tsv --output-fp /tmpfs/tmp.biom --sc-separated "taxonomy"
    biom convert --input-fp /tmpfs/tmp.biom --output-fp $3 --to-tsv --process-obs-metadata "taxonomy" --tsv-metadata-formatter "sc_separated" --header-key "taxonomy"
else
    echo "This script will work if and only if it has **three** arguments"
fi
