# This python script is developed by Chunchao Zhang and Bicna Song for CNH IMPACT trial

import os
import csv
import sys
import re
import math
import itertools
import argparse
import gzip
from pathlib import Path
from util.digestion import tryptic_pep, lysC_pep

def read_tsv(file_path):
    with open(file_path, "r") as f:
        data = csv.reader(f, delimiter="\t")
        allMut = [row for row in data]
    return allMut

def load_proteomes():
    root_dir = next(
        p for p in [Path.cwd()] + list(Path.cwd().parents) if (p / ".git").is_dir()
    )
    data_dir = root_dir / "data"

    ensembl = data_dir / "Homo_sapiens.GRCh38.pep.all.fa.gz"
    refseq = data_dir / "GCF_000001405.40_GRCh38.p14_protein.faa.gz"

    ensp = re.compile(r"(ENSP[0-9a-zA-Z]+)")
    ensp2pep = {}
    with gzip.open(ensembl, "rt") as f:
        for line in f:
            if line.startswith(">"):
                ENSP_id = ensp.findall(line)[0]
                Seq = ""
            else:
                Seq += line.rstrip("\n")
                ensp2pep[ENSP_id] = Seq

    refseqID2prot = {}
    with gzip.open(refseq, "rt") as f:
        for line in f:
            if line.startswith(">"):
                refseq_id = line.split()[0][1:].split(".")[0]
                Seq = ""
            else:
                Seq += line.rstrip("\n")
                refseqID2prot[refseq_id] = Seq

    return ensp2pep, refseqID2prot

def mutant_pept(file_path, digestion_func, ensp2pep, refseqID2prot):
    Mutations = read_tsv(file_path)
    snv = [
        row for row in Mutations
        if row[8] in {"Missense_Mutation", "Nonsense_Mutation"} and row[9] == "SNP"
    ]
    header = Mutations[0]

    tumor_idx = header.index("Tumor_Sample_Barcode")
    norm_idx = header.index("Matched_Norm_Sample_Barcode")
    GN_idx = header.index("Hugo_Symbol")
    gene_idx = header.index("Gene")
    trans_idx = header.index("Transcript_ID")
    prot_idx = header.index("ENSP")
    genVar_idx = header.index("HGVSc")
    protVar_idx = header.index("HGVSp_Short")
    mutType_idx = header.index("Variant_Classification")
    VarType_idx = header.index("Variant_Type")

    output = []
    AA_wt = re.compile(r"p\.([a-zA-Z]{1})\d")
    Mut_pos = re.compile(r"\d+")
    AA_mut = re.compile(r"\d+(.*)")
    i = 0

    for line in snv:
        Barcode = line[tumor_idx] + ":" + line[norm_idx]
        GN = line[GN_idx]
        geneID = line[gene_idx]
        transcriptID = line[trans_idx]
        protID = line[prot_idx]
        genVar = line[genVar_idx]
        protVar = line[protVar_idx]
        mutType = line[mutType_idx]
        VarType = line[VarType_idx]

        try:
            get_wt = re.findall(AA_wt, protVar)[0]
            get_mut = re.findall(AA_mut, protVar)[0]
            get_pos = int(re.findall(Mut_pos, protVar)[0])
        except:
            continue

        if protID.startswith("ENSP"):
            if protID not in ensp2pep:
                continue
            protein = ensp2pep[protID]
        else:
            new_protID = protID.split(".")[0]
            if new_protID not in refseqID2prot:
                continue
            protein = refseqID2prot[new_protID]

        if get_pos > len(protein):
            continue
        if protein[get_pos - 1] != get_wt:
            continue

        mut_protein = protein[0:get_pos - 1] + get_mut + protein[get_pos:]
        mut_pep = digestion_func(get_pos, mut_protein)
        i += 1
        peptID = "snv_" + str(i)
        mut_pep_id = (
            ">" + peptID + "|" + Barcode + "|" + GN + "|" + geneID + "|" +
            transcriptID + "|" + protID + "|" + genVar + "|" + protVar + "|" +
            mutType + "|peptMutPos=" + str(mut_pep[0])
        )
        output.append([mut_pep_id, mut_pep[1].split("*")[0]])

    return output

def write_fasta(input_file, output_file, digestion_func, ensp2pep, refseqID2prot):
    mut_pept = mutant_pept(input_file, digestion_func, ensp2pep, refseqID2prot)
    with open(output_file, "w") as f:
        for i in mut_pept:
            f.write(i[0] + "\n")
            f.write(i[1] + "\n")

def main():
    parser = argparse.ArgumentParser(description="A python script for SNVs")
    parser.add_argument("-i", "--input", required=True, help="Path to input file")
    parser.add_argument("-o", "--output", required=True, help="Path to output file")
    parser.add_argument(
        "-d", "--digestion",
        choices=["trypsin", "lysC"],
        required=True,
        help="Enzymes: trypsin or lysC"
    )
    if "--help" in sys.argv or "-h" in sys.argv:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    ensp2pep, refseqID2prot = load_proteomes()
    digestion_func = tryptic_pep if args.digestion == "trypsin" else lysC_pep
    write_fasta(args.input, args.output, digestion_func, ensp2pep, refseqID2prot)
    print(f"✅ Done: {args.output}")

if __name__ == "__main__":
    main()
