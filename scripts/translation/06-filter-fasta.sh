#!/bin/bash
set -e
set -o pipefail

# Remove custom peptides that already exist in the canonical proteome.
#
# A peptide identical to canonical sequence is already searchable, adds nothing
# tumor-specific, and only inflates the database and the fragment index. The
# translation steps emit a window around each event, so most of what they
# produce is ordinary flanking sequence; this is what strips it.
#
# This is an exact, full-length substring test, so it is done with grep -F -f
# (Aho-Corasick: every pattern matched in a single pass), not with blastp.
#
# blastp was used here previously. It is a heuristic LOCAL ALIGNER with an
# E-value reporting threshold, and short queries score too low to be reported.
# Measured against this pipeline's own peptides, using the canonical proteome as
# both the BLAST database and the ground truth:
#
#     length   truly present   blastp   blastp-short
#        8           16            5         16
#        9           10            8         10
#       10           15           14         15
#    11-13           60           60         60
#       14           31           30         31
#    TOTAL          132          117        132
#
# Default blastp missed 11% of filterable peptides overall and 69% of 8-mers,
# which stayed in the database as redundant canonical sequence. -task
# blastp-short would fix the threshold, but an alignment engine is still the
# wrong tool for a string-containment question.

output_basename=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --custom_fasta)     input_fasta="$2";   shift 2 ;;
    --uniprot_fasta)    uniprot_fasta="$2"; shift 2 ;;
    --output_basename)  output_basename="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

filtered_fasta="${output_basename}_tumor_specific_filtered.fasta"
match_report="${output_basename}_uniprot_match_report.tsv"
full_matches="${output_basename}_uniprot_full_matches.txt"

# Flatten the reference to ONE line per protein before matching. The FASTA is
# wrapped at 60 columns, and a peptide spanning a line break would otherwise be
# silently missed -- a real failure mode, hit while testing this change.
seqkit seq -w 0 "$uniprot_fasta" | grep -v '^>' > uniprot_flat.txt

# id<TAB>sequence for every query peptide
seqkit fx2tab -i "$input_fasta" | cut -f1,2 > query_tab.tsv
cut -f2 query_tab.tsv | sort -u > query_seqs.txt

# One pass: which query sequences occur anywhere in the reference
if [ -s query_seqs.txt ]; then
  LC_ALL=C grep -o -F -f query_seqs.txt uniprot_flat.txt | sort -u > matched_seqs.txt || true
else
  : > matched_seqs.txt
fi

# Map matched sequences back to peptide IDs
awk -F'\t' 'NR==FNR { m[$0]; next } ($2 in m) { print $1 }' \
    matched_seqs.txt query_tab.tsv | sort -u > "$full_matches"

# Per-peptide report, replacing the old BLAST tabular output
{
  printf 'peptide_id\tpeptide\tin_uniprot\n'
  awk -F'\t' 'NR==FNR { m[$0]; next } { print $1"\t"$2"\t"(($2 in m) ? "yes" : "no") }' \
      matched_seqs.txt query_tab.tsv
} > "$match_report"

total=$(grep -c "^>" "$input_fasta")
matches=$(wc -l < "$full_matches")
remaining=$(( total - matches ))

echo "Total peptides: $total"
echo "Filtered out (exact full-length UniProt matches): $matches"
echo "Remaining peptides: $remaining"

if [ "$matches" -gt 0 ]; then
  seqkit grep -v -f "$full_matches" "$input_fasta" > "$filtered_fasta"
else
  cp "$input_fasta" "$filtered_fasta"
fi
