#!/bin/bash

# Author: Bicna Song
# Description: Run translation scripts for SNV and fusion inputs

set -e
set -o pipefail

# Set up paths
root_dir=$(git rev-parse --show-toplevel)
data_dir="$root_dir/data"
analysis_dir="$root_dir/analyses/translation"
results_dir="$analysis_dir/results"

# Digestion enzyme (choose: lysC or trypsin)
digestion="lysC"

# Create results directory if it does not exist
mkdir -p "$results_dir"

# Run SNV translation
echo "Running SNV translation with digestion..."
python "$analysis_dir/01-snv-translate.py" \
  -i "$data_dir/snv-consensus-plus-hotspots.maf-test.tsv" \
  -o "$results_dir/snv-consensus-plus-hotspots.maf.fasta" \
  -d "$digestion"

# Run STAR fusion translation from gzipped input using pipe
echo "Running fusion translation with digestion from compressed input.."
gunzip -c "$data_dir/fusion-starfusion-test.tsv.gz" | \
  python "$analysis_dir/02-star-fusion-translate.py" \
    -i - \
    -o "$results_dir/fusion-starfusion-peptides.fasta" \
    -d "$digestion"

# Run Arriba fusion translation from gzipped input using pipe
echo "Running fusion translation with digestion from compressed input.."
gunzip -c "$data_dir/fusion-arriba-test.tsv.gz" | \
  python "$analysis_dir/03-arriba-fusion-translate.py" \
    -i - \
    -o "$results_dir/fusion-arriba-peptides.fasta" \
    -d "$digestion"

echo "✅ Translation completed. Outputs saved to $results_dir"
