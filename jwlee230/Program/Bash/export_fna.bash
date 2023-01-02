#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

if (($# == 2)); then
    qiime tools export --input-path $1 --output-path /tmpfs
    mv /tmpfs/dna-sequences.fasta $2
else
    echo "This script will work if and only if it has **two** arguments"
fi
