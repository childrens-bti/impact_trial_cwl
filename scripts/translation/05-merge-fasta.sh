#!/bin/bash

set -e
set -o pipefail

# Define the output file
output_file="results/custom.fasta"

# Create an empty file or clear the existing file
> "$output_file"

# Loop through all *peptides.fasta files directly in results/
for fasta_file in results/*peptides.fasta; do
    if [ -f "$fasta_file" ]; then
        cat "$fasta_file" >> "$output_file"
    fi
done
