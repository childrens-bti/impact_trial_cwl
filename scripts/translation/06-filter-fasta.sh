#!/bin/bash
set -e
set -o pipefail

# Parse command line arguments
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
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

filtered_fasta=custom_filtered.fasta

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
       -out blast_results.tsv \
       -outfmt "6 qseqid sseqid pident length qlen slen"

# extract query IDs with full-length 100% matches
awk '$3 == 100 && $4 == $5 {print $1}' blast_results.tsv | sort -u > full_matches.txt

# count matches and total peptides
total=$(grep -c "^>" "$input_fasta")
matches=$(wc -l < full_matches.txt)
remaining=$(( total - matches ))

echo "Total peptides: $total"
echo "Filtered out (100% UniProt matches): $matches"
echo "Remaining peptides: $remaining"

# filter out 100% match peptides 
seqkit grep -v -f full_matches.txt "$input_fasta" \
    > "$filtered_fasta"
