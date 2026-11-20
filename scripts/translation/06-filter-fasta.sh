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
# This is an exact, full-length substring test: each unique peptide is checked
# on its own for containment in the reference (Python `in`), not with blastp.
#
# Do NOT replace this with `grep -o -F -f peptides ref`. grep -o reports only
# non-overlapping matches (leftmost-longest at each position, then it resumes
# after the match), so a peptide that overlaps or is a prefix of another
# matched peptide is never reported. The translation steps emit overlapping
# windows, so this would leave canonical peptides unfiltered.
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

# Per-peptide containment check. Both FASTAs are parsed record by record, so
# line wrapping in the reference (60 columns) cannot hide a peptide that spans
# a line break. Proteins are joined with a newline, which never occurs in a
# peptide, so a match cannot straddle two proteins.
python3 - "$input_fasta" "$uniprot_fasta" "$match_report" "$full_matches" <<'PY'
import sys

query_fasta, ref_fasta, report_path, matches_path = sys.argv[1:5]

def read_fasta(path):
    """Yield (id, sequence); id is the header up to the first whitespace,
    matching how seqkit identifies records."""
    rid, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(chunks)
                rid, chunks = (line[1:].split() or [""])[0], []
            else:
                chunks.append(line)
    if rid is not None:
        yield rid, "".join(chunks)

reference = "\n".join(seq for _, seq in read_fasta(ref_fasta))
queries = list(read_fasta(query_fasta))

# An empty sequence is trivially "in" any string; never treat it as a match.
in_ref = {seq: bool(seq) and seq in reference for seq in {seq for _, seq in queries}}

with open(report_path, "w") as report:
    report.write("peptide_id\tpeptide\tin_uniprot\n")
    for pid, seq in queries:
        report.write(f"{pid}\t{seq}\t{'yes' if in_ref[seq] else 'no'}\n")

matched_ids = sorted({pid for pid, seq in queries if in_ref[seq]})
with open(matches_path, "w") as out:
    out.writelines(f"{pid}\n" for pid in matched_ids)
PY

# grep -c exits 1 when the count is 0; do not let set -e kill an empty input
total=$(grep -c "^>" "$input_fasta" || true)
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
