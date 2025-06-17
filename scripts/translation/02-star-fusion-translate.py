# This python script is developed by Chunchao Zhang and Bicna Song for CNH IMPACT trial

import os
import csv
import sys
import re
import math
import argparse
from util.digestion import tryptic_pep, lysC_pep

def fusion_prot(file_path, digestion_func):
    """Outputs fusion protein sequences and their FASTA headers based on fusion type."""
    allFusion = read_tsv(file_path)

    # star fusion table header
    header = allFusion[0]
    
    # get column numbers for different variables
    idx = {
        col: header.index(col)
        for col in [
            "tumor_id",
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
    count = 0

    # Loop through all fusion events
    for line in allFusion[1:]:
       fusion_type = line[idx["PROT_FUSION_TYPE"]]
       if fusion_type not in {"INFRAME", "FRAMESHIFT"}:
            continue
       
       count += 1
       full_seq = line[idx["FUSION_TRANSL"]].split("*")[0]

        # Determine fusion site (only relevant for INFRAME)
       try:
           fusion_site = math.ceil(int(line[idx["CDS_LEFT_RANGE"]].split("-")[1]) / 3)
       except:
            fusion_site = None

        # Apply digestion for INFRAME fusions
       if fusion_type == "INFRAME" and fusion_site is not None:
            try:
               _, peptide = digestion_func(fusion_site, full_seq)
            except Exception:
                peptide = full_seq  # fallback
       else:
           peptide = full_seq

       # Compose header
       fusion_seq_id = (
           f">star_fusion_{count}|{line[idx['tumor_id']]}|{line[idx['FusionName']]}|"
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

def write_fasta(input_file, output_file, digestion_func):
    """output fasta file"""
    fusion_seq = fusion_prot(input_file, digestion_func)
    with open(output_file, "w") as f:
        for i in fusion_seq:
            f.write(i[0] + "\n")
            f.write(i[1] + "\n")

def main():
    """Main function with arguments"""
    parser = argparse.ArgumentParser(description="A python script for star fusions")
    # Adding arguments
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to fusion TSV file")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to output FASTA file")
    parser.add_argument("-d", "--digestion", choices=["trypsin", "lysC"], default="lysC",
                        help="Protease used for digestion (default: lysC)")
    # Parsing arguments
    args = parser.parse_args()

    digestion_func = tryptic_pep if args.digestion == "trypsin" else lysC_pep

    # Process the file and write the results
    write_fasta(args.input, args.output, digestion_func)
    # Verbosity
    print(f"Processing data from {args.input} and saving results to {args.output}")

if __name__ == "__main__":
    main()
