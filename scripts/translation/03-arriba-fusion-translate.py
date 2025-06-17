# This python script is developed by Chunchao Zhang and Bicna Song for CNH IMPACT trial

import os
import csv
import sys
import re
import math
import argparse
from util.digestion import tryptic_pep, lysC_pep

def clean_peptide(pep):
    """
    Convert to uppercase and remove any character that is not one of the 20 standard amino acid codes.
    """
    pep_upper = pep.upper()
    cleaned = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", pep_upper)
    return cleaned

def fusion_prot(file_path, digestion_func):
    """
    Outputs fusion protein sequences and their FASTA headers based on fusion type.
    """
    allFusion = read_tsv(file_path)

    # arriba fusion table header
    header = allFusion[0]
    
    # get column numbers for different variables
    idx = {
        col: header.index(col)
        for col in [
            "tumor_id",
            "gene1--gene2",
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
    count = 0

    # Loop through all fusion events
    for line in allFusion[1:]:
       fusion_type = line[idx["reading_frame"]]
       if fusion_type not in {"in-frame", "out-of-frame"}:
            continue
       
       seq = line[idx["peptide_sequence"]]

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

       # Apply digestion for fusions
       new_pos, pep_cut = digestion_func(fusion_pos, full_seq)

       # Compose header
       count += 1
       fusion_seq_id = (
           f">arriba_fusion_{count}|{line[idx['tumor_id']]}|{line[idx['gene1--gene2']]}|"
           f"{line[idx['gene1']]}|{line[idx['breakpoint1']]}|{line[idx['transcript_id1']]}|"
           f"{line[idx['gene2']]}|{line[idx['breakpoint2']]}|{line[idx['transcript_id2']]}|"
           f"{fusion_type}|Protein_fusion_site:{fusion_pos}"
        )
       
       output.append([fusion_seq_id, pep_cut])
       
    return output

def read_tsv(file_path):
    if file_path == "-":
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
