# This python script is developed by Chunchao Zhang and Bicna Song for CNH IMPACT trial

import os
import csv
import sys
import re
import math
import argparse
import gzip
from util.digestion import tryptic_pep, lysC_pep
from util.annofuse_filtering import load_annofuse, filter_starfusion_by_annofuse

def normalize_star_header(rows):
    """
    Ensure STAR-Fusion header uses 'FusionName' (not '#FusionName').
    """
    if not rows:
        return rows
    header = rows[0][:]
    for i, col in enumerate(header):
        if col.strip() == "#FusionName":
            header[i] = "FusionName"
    rows[0] = header
    return rows

def load_gtex_pairs(path):
    """Load ordered FusionName strings from GTEx TSV (must have 'FusionName')"""
    normals = set()
    with gzip.open(path, "rt", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "FusionName" not in reader.fieldnames:
            raise ValueError("GTEx normals must contain column 'FusionName'")
        for row in reader:
            val = (row.get("FusionName") or "").strip()
            if val:
                normals.add(val.upper())
    print(f"Loaded {len(normals)} GTEx normal fusion entries from {path}")
    return normals
  
def filter_star_by_gtex_normals(allFusion_rows, normals_set):
    """Keep header + rows whose FusionName are NOT in GTEx normals."""
    if not allFusion_rows or not normals_set:
        return allFusion_rows
    header = allFusion_rows[0]
    iF = header.index("FusionName")
    kept, removed = [header], 0
    for row in allFusion_rows[1:]:
        name = (row[iF] or "").strip().upper()
        if name in normals_set:
            removed += 1
            continue
        kept.append(row)
    print(f"GTEx ordered filter: removed {removed}, kept {len(kept)-1}")
    return kept

def fusion_prot(file_path, digestion_func, skip_digestion, gtex_normals=None, annofuse_filter=None):
    """Outputs fusion protein sequences and their FASTA headers based on fusion type."""
    allFusion = read_tsv(file_path)
    allFusion = normalize_star_header(allFusion)
    
    # Drop GTEx-normal fusions
    if gtex_normals:
      normals = load_gtex_pairs(gtex_normals)
      allFusion = filter_star_by_gtex_normals(allFusion, normals)

    # Keep only patient AnnoFuse calls
    if annofuse_filter:
        keep_set = load_annofuse(annofuse_filter)
        allFusion = filter_starfusion_by_annofuse(allFusion, keep_set)

    # star fusion table header
    header = allFusion[0]
    
    # get column numbers for different variables
    idx = {
        col: header.index(col)
        for col in [
            "FusionName",
            "LeftGene",
            "LeftBreakpoint",
            "RightGene",
            "RightBreakpoint",
            "CDS_LEFT_ID",
            "CDS_LEFT_RANGE",
            "CDS_RIGHT_ID",
            "CDS_RIGHT_RANGE",
            "PROT_FUSION_TYPE",
            "FUSION_TRANSL",
        ]
    }

    output = []

    # Loop through all fusion events
    for line in allFusion[1:]:
       fusion_type = line[idx["PROT_FUSION_TYPE"]]
       if fusion_type not in {"INFRAME", "FRAMESHIFT"}:
            continue
       
       full_seq = line[idx["FUSION_TRANSL"]].split("*")[0]

        # Determine fusion site (only relevant for INFRAME)
       try:
           fusion_site = math.ceil(int(line[idx["CDS_LEFT_RANGE"]].split("-")[1]) / 3)
       except:
            fusion_site = None

        # Apply digestion only if not skipping and fusion is in-frame
       if not skip_digestion and fusion_type == "INFRAME" and fusion_site is not None:
            try:
                _, peptide = digestion_func(fusion_site, full_seq)
            except Exception:
                peptide = full_seq  # fallback on error
       else:
            peptide = full_seq

       # Compose header
       fusion_seq_id = (
           f">star_fusion|{line[idx['FusionName']]}|"
           f"{line[idx['LeftGene']]}|{line[idx['LeftBreakpoint']]}|{line[idx['RightGene']]}|"
           f"{line[idx['RightBreakpoint']]}|leftCDS_id:{line[idx['CDS_LEFT_ID']]}|{line[idx['CDS_LEFT_RANGE']]}|"
           f"rightCDS_id:{line[idx['CDS_RIGHT_ID']]}|{line[idx['CDS_RIGHT_RANGE']]}|"
           f"{fusion_type}|Protein_fusion_site:{fusion_site}"
        )
       
       output.append([fusion_seq_id, peptide])
       
    return output

def read_tsv(file_path):
    if file_path == "-":  # Support stdin
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
    parser = argparse.ArgumentParser(description="A python script for STAR fusions")
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
                        help="Path to GTEx normals TSV with column 'FusionName'")
    parser.add_argument("--annofuse_filter", type=str, default=None,
                        help="Path to *.annoFuse_filter.tsv (must contain 'FusionName'); "
                             "only these fusions will be translated")

    # Parsing arguments
    args = parser.parse_args()

    digestion_func = tryptic_pep if args.digestion == "trypsin" else lysC_pep

    # Generate the list and always materialize output for CWL collection
    records = fusion_prot(args.input, digestion_func, args.skip_digestion, args.gtex_normals, args.annofuse_filter)

    # Always write the FASTA, even when empty
    write_fasta(args.output, records)

    if not records:
        print("No INFRAME or FRAMESHIFT entries found; wrote empty FASTA.")
        sys.exit(0)

    print(f"Processing data from {args.input} and saving results to {args.output}")
    if args.skip_digestion:
        print("Note: digestion step was skipped.")
    if args.gtex_normals:
        print(f"Filtered by GTEx STAR-Fusion output: {args.gtex_normals}")
    if args.annofuse_filter:
        print(f"Filtered by Patient AnnoFuse output: {args.annofuse_filter}")

if __name__ == "__main__":
    main()
