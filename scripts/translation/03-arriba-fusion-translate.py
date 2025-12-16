# This python script is developed by Chunchao Zhang and Bicna Song for CNH IMPACT trial

import os
import csv
import sys
import re
import math
import argparse
import gzip
from util.digestion import tryptic_pep, lysC_pep
from util.annofuse_filtering import load_annofuse, filter_arriba_by_annofuse

def clean_peptide(pep):
    """
    Convert to uppercase and remove any character that is not one of the 20 standard amino acid codes.
    """
    pep_upper = pep.upper()
    cleaned = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", pep_upper)
    return cleaned
  
def normalize_arriba_header(rows):
    """
    Ensure Arriba header uses 'gene1' (not '#gene1').
    Works on list-of-lists returned by read_tsv.
    """
    if not rows:
        return rows
    header = rows[0][:]
    for i, col in enumerate(header):
        if col.strip() == "#gene1":
            header[i] = "gene1"
    rows[0] = header
    return rows

def load_gtex_pairs(path):
    """Load ordered 'GENE1--GENE2' keys from GTEx TSV (mush have 'gene1--gene2')"""
    normals = set()
    with gzip.open(path, "rt", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "gene1--gene2" not in reader.fieldnames:
            raise ValueError("GTEx normals must contain column 'gene1--gene2'")
        for row in reader:
            val = (row.get("gene1--gene2") or "").strip()
            if val:
                normals.add(val.upper())
    print(f"Loaded {len(normals)} GTEx normal fusions from {path}")
    return normals

def filter_arriba_by_gtex_normals(allFusion_rows, normals_set):
    """Keep header + rows whose ordered (gene1,gene2) are NOT in GTEx normals."""
    if not allFusion_rows or not normals_set:
        return allFusion_rows
    header = allFusion_rows[0]
    i1 = header.index("gene1")
    i2 = header.index("gene2")
    kept, removed = [header], 0
    for row in allFusion_rows[1:]:
        key = f"{row[i1].strip().upper()}--{row[i2].strip().upper()}"
        if key in normals_set:
            removed += 1
            continue
        kept.append(row)
    print(f"GTEx ordered filter: removed {removed}, kept {len(kept)-1}")
    return kept

def fusion_prot(file_path, digestion_func, skip_digestion, gtex_normals=None, annofuse_filter=None):
    """
    Outputs fusion protein sequences and their FASTA headers based on fusion type.
    """
    allFusion = read_tsv(file_path)
    allFusion = normalize_arriba_header(allFusion)
    
    # Drop GTEx-normal fusions
    if gtex_normals:
      normals = load_gtex_pairs(gtex_normals)
      allFusion = filter_arriba_by_gtex_normals(allFusion, normals)

    # Keep only patient AnnoFuse calls
    if annofuse_filter:
        keep_set = load_annofuse(annofuse_filter)
        allFusion = filter_arriba_by_annofuse(allFusion, keep_set)

    # arriba fusion table header
    header = allFusion[0]
    
    # get column numbers for different variables
    idx = {
        col: header.index(col)
        for col in [
            "gene1",
            "breakpoint1",
            "gene2",
            "breakpoint2",
            "transcript_id1",
            "transcript_id2",
            "reading_frame",
            "peptide_sequence",
        ]
    }

    output = []

    # Loop through all fusion events
    for line in allFusion[1:]:
       fusion_type = line[idx["reading_frame"]]
       if fusion_type not in {"in-frame", "out-of-frame"}:
            continue
       
       seq = line[idx["peptide_sequence"]].split("*")[0]

       if "|" in seq:
            parts = seq.split("|", 1)
            raw_left, raw_right = parts[0], parts[1]
       else:
            # No pipe found: treat entire seq as left side, right side = empty
            raw_left = seq
            raw_right = ""

       left_clean = clean_peptide(raw_left)
       right_clean = clean_peptide(raw_right)

       fusion_pos = len(left_clean)
       full_seq = left_clean + right_clean

       if not skip_digestion:
           try:
               _, peptide = digestion_func(fusion_pos, full_seq)
           except Exception:
               peptide = full_seq
       else:
           peptide = full_seq

       # Compose header
       fusion_name = f"{line[idx['gene1']]}--{line[idx['gene2']]}"
       fusion_seq_id = (
           f">arriba_fusion|{fusion_name}|"
           f"{line[idx['gene1']]}|{line[idx['breakpoint1']]}|{line[idx['transcript_id1']]}|"
           f"{line[idx['gene2']]}|{line[idx['breakpoint2']]}|{line[idx['transcript_id2']]}|"
           f"{fusion_type}|Protein_fusion_site:{fusion_pos}"
        )
       output.append([fusion_seq_id, peptide])
       
    return output

def read_tsv(file_path):
    if file_path == "-":
        f = sys.stdin
    else:
        f = open(file_path, 'r')
    with f:
        data = csv.reader(f, delimiter='\t')
        return [row for row in data]

def write_fasta(output_file, records):
    """Output fasta file"""
    with open(output_file, "w") as f:
        for header, seq in records:
            f.write(header + "\n")
            f.write(seq + "\n")

def main():
    """Main function with arguments"""
    parser = argparse.ArgumentParser(description="A python script for Arriba fusions")
    # Adding arguments
    parser.add_argument("-i", "--input", type=str, required=True, 
                        help="Path to fusion TSV file")
    parser.add_argument("-o", "--output", type=str, required=True, 
                        help="Path to output FASTA file")
    parser.add_argument("-d", "--digestion", choices=["trypsin", "lysC"], default="lysC",
                        help="Protease used for digestion (default: lysC)")
    parser.add_argument("--skip_digestion", action="store_true",
                        help="If set, skip the peptide digestion step entirely")
    parser.add_argument("--gtex_normals", type=str, default=None,
                        help="Path to GTEx normals TSV with an ordered column 'gene1--gene2'")
    parser.add_argument("--annofuse_filter", type=str, default=None,
                        help="Path to *.annoFuse_filter.tsv (must contain 'FusionName'); "
                             "only these fusions will be translated")

    # Parsing arguments
    args = parser.parse_args()

    digestion_func = tryptic_pep if args.digestion == "trypsin" else lysC_pep

    records = fusion_prot(args.input, digestion_func, args.skip_digestion, args.gtex_normals, args.annofuse_filter)
    
    if not records:
        print("No in-frame or out-of-frame entries found; nothing to translate.")
        sys.exit(0)

    write_fasta(args.output, records)
    print(f"Processing data from {args.input} and saving results to {args.output}")
    if args.skip_digestion:
        print("Note: digestion step was skipped.")
    if args.gtex_normals:
        print(f"Filtered by GTEx Arriba output: {args.gtex_normals}")
    if args.annofuse_filter:
        print(f"Filtered by Patient AnnoFuse output: {args.annofuse_filter}")

if __name__ == "__main__":
    main()
