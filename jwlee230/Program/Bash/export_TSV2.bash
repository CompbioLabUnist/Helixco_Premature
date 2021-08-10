#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    qiime tools export --input-path $1 --output-path /tmpfs
    biom convert --input-fp /tmpfs/feature-table.biom --output-fp $2 --to-tsv
    sed --in-place '1d' $2
else
    echo "This script will work if and only if it has **two** arguments"
fi

