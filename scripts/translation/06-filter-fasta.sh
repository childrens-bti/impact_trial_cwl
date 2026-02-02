#!/bin/bash
set -e
set -o pipefail

# Parse command line arguments
dna_sample_id=""
rna_sample_id=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --custom_fasta)
      input_fasta="$2"
      shift 2
      ;;
    --uniprot_fasta)
      uniprot_fasta="$2"
      shift 2
      ;;
    --dna_sample_id)
      dna_sample_id="$2"
      shift 2
      ;;
    --rna_sample_id)
      rna_sample_id="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Generate output filenames with sample IDs
filtered_fasta="${dna_sample_id}_${rna_sample_id}_tumor_specific_filtered.fasta"
blast_results="${dna_sample_id}_${rna_sample_id}_uniprot_blast_results.tsv"
full_matches="${dna_sample_id}_${rna_sample_id}_uniprot_full_matches.txt"

# if BLAST database not found, build
if [ ! -f "uniprot_db.pin" ]; then
    echo "BLAST database not found. Building now..."
    makeblastdb -in "$uniprot_fasta" -dbtype prot -out uniprot_db
else
    echo "BLAST database already exists. Skipping build."
fi

# run blastp
blastp -query "$input_fasta" \
       -db uniprot_db \
       -out "$blast_results" \
       -outfmt "6 qseqid sseqid pident length qlen slen"

# extract query IDs with full-length 100% matches
awk '$3 == 100 && $4 == $5 {print $1}' "$blast_results" | sort -u > "$full_matches"

# count matches and total peptides
total=$(grep -c "^>" "$input_fasta")
matches=$(wc -l < "$full_matches")
remaining=$(( total - matches ))

echo "Total peptides: $total"
echo "Filtered out (100% UniProt matches): $matches"
echo "Remaining peptides: $remaining"

# filter out 100% match peptides 
seqkit grep -v -f "$full_matches" "$input_fasta" \
    > "$filtered_fasta"
