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
from Bio import SeqIO
from Bio.Seq import Seq
from util.digestion import tryptic_pep, lysC_pep


def read_tsv(file_path):
    """
    Read a TSV file, skipping any lines starting with '#' (e.g., version lines).
    """
    with open(file_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        # drop comment lines (first cell starts with '#')
        return [row for row in reader if row and not row[0].startswith("#")]


def load_proteomes():
    root_dir = next(
        p for p in [Path.cwd()] + list(Path.cwd().parents) if (p / ".git").is_dir()
    )
    data_dir = root_dir / "data"
    ensembl = data_dir / "Homo_sapiens.GRCh38.pep.all.fa.gz"
    refseq = data_dir / "GCF_000001405.40_GRCh38.p14_protein.faa.gz"
    cds_file = data_dir / "Homo_sapiens.GRCh38.cds.all.fa.gz"

    ensp2pep = {}
    refseqID2prot = {}
    transcript2cds = {}

    with gzip.open(ensembl, "rt") as f:
        for line in f:
            if line.startswith(">"):
                ensp_id = re.findall(r"(ENSP[0-9a-zA-Z]+)", line)[0]
                Seq = ""
            else:
                Seq += line.strip()
                ensp2pep[ensp_id] = Seq

    with gzip.open(refseq, "rt") as f:
        for line in f:
            if line.startswith(">"):
                refseq_id = line.split()[0][1:].split(".")[0]
                Seq = ""
            else:
                Seq += line.strip()
                refseqID2prot[refseq_id] = Seq

    with gzip.open(cds_file, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            transcript_id = record.id.split("|")[0].split(".")[0]
            transcript2cds[transcript_id] = str(record.seq)

    return ensp2pep, refseqID2prot, transcript2cds


def parse_inframe_variant(protVar):
    match_delins = re.match(r"p\.([A-Z])(\d+)_([A-Z])(\d+)delins([A-Z]+)", protVar)
    if match_delins:
        return {
            "type": "delins",
            "start": int(match_delins.group(2)),
            "end": int(match_delins.group(4)),
            "inserted": match_delins.group(5),
        }

    match_del = re.match(r"p\.([A-Z])(\d+)del", protVar)
    if match_del:
        pos = int(match_del.group(2))
        return {"type": "del", "start": pos, "end": pos, "inserted": ""}

    match_del_multi = re.match(r"p\.([A-Z])(\d+)_([A-Z])(\d+)del", protVar)
    if match_del_multi:
        return {
            "type": "del",
            "start": int(match_del_multi.group(2)),
            "end": int(match_del_multi.group(4)),
            "inserted": "",
        }

    match_dup = re.match(r"p\.([A-Z])(\d+)dup", protVar)
    if match_dup:
        pos = int(match_dup.group(2))
        return {"type": "dup", "start": pos, "end": pos, "inserted": ""}

    match_dup_multi = re.match(r"p\.([A-Z])(\d+)_([A-Z])(\d+)dup", protVar)
    if match_dup_multi:
        return {
            "type": "dup",
            "start": int(match_dup_multi.group(2)),
            "end": int(match_dup_multi.group(4)),
            "inserted": "",
        }

    match_ins = re.match(r"p\.([A-Z])(\d+)_([A-Z])(\d+)ins([A-Z]+)", protVar)
    if match_ins:
        return {
            "type": "ins",
            "start": int(match_ins.group(2)),
            "end": int(match_ins.group(4)),
            "inserted": match_ins.group(5),
        }

    return None


def parse_hgvsc_variant(hgvsc):
    # Handle insertions
    if "ins" in hgvsc:
        match = re.match(r"c\.(\d+)_?[0-9+]*ins([ACGT]+)", hgvsc)
        if match:
            start = int(match.group(1))
            seq = match.group(2)
            return ("ins", start, seq)

    # Handle deletions
    elif "del" in hgvsc:
        match = re.match(r"c\.(\d+)(?:_(\d+))?del", hgvsc)
        if match:
            start = int(match.group(1))
            end = int(match.group(2)) if match.group(2) else start
            return ("del", start, end)

    # Handle duplications
    elif "dup" in hgvsc:
        match = re.match(r"c\.(\d+)(?:_(\d+))?dup", hgvsc)
        if match:
            start = int(match.group(1))
            end = int(match.group(2)) if match.group(2) else start
            return ("dup", start, end)

    return None


def translate_frameshift_cds(cds_seq, var_type, start, info):
    """
    Translate mutant CDS resulting from frameshift-inducing insertion, deletion, or duplication.
    Truncates protein at first stop codon (*).

    Parameters:
        cds_seq (str): Wild-type CDS sequence
        var_type (str): One of "ins", "del", or "dup"
        start (int): 1-based coordinate for mutation
        info (str or int): inserted bases (if ins), or end position (if del/dup)

    Returns:
        str: Mutant protein sequence, truncated at first stop codon
    """
    if var_type == "ins":
        mutant_cds = cds_seq[:start] + info + cds_seq[start:]
    elif var_type == "del":
        end = info
        mutant_cds = cds_seq[: start - 1] + cds_seq[end:]
    elif var_type == "dup":
        end = info
        dup_seq = cds_seq[start - 1 : end]
        mutant_cds = cds_seq[:end] + dup_seq + cds_seq[end:]
    else:
        return None

    # Trim if sequence length is not divisible by 3
    trim_len = len(mutant_cds) % 3
    if trim_len != 0:
        mutant_cds = mutant_cds[: len(mutant_cds) - trim_len]

    # Translate
    protein = str(Seq(mutant_cds).translate(to_stop=False))
    return protein


def mutant_pept(
    file_path, digestion_func, ensp2pep, refseqID2prot, transcript2cds, skip_digestion
):
    Mutations = read_tsv(file_path)
    if not Mutations:
        return []

    header = Mutations[0]
    idx = {col: i for i, col in enumerate(header)}

    valid_snv = {"Missense_Mutation", "Nonsense_Mutation"}
    valid_indel = {"In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins"}

    output = []

    for row in Mutations[1:]:
        mutType = row[idx["Variant_Classification"]]
        varType = row[idx["Variant_Type"]]
        if (mutType in valid_snv and varType == "SNP") or (
            mutType in valid_indel and varType in {"DEL", "INS"}
        ):
            Barcode = (
                row[idx["Tumor_Sample_Barcode"]]
                + ":"
                + row[idx["Matched_Norm_Sample_Barcode"]]
            )
            GN = row[idx["Hugo_Symbol"]]
            geneID = row[idx["Gene"]]
            transcriptID = row[idx["Transcript_ID"]].split(".")[0]
            protID = row[idx["ENSP"]]
            genVar = row[idx["HGVSc"]]
            protVar = row[idx["HGVSp_Short"]]

            # Load protein sequence
            if protID.startswith("ENSP"):
                if protID not in ensp2pep:
                    continue
                protein = ensp2pep[protID]
            else:
                ref_id = protID.split(".")[0]
                if ref_id not in refseqID2prot:
                    continue
                protein = refseqID2prot[ref_id]

            # SNVs
            if mutType in valid_snv:
                match = re.match(r"p\.([A-Z])(\d+)([A-Z*])", protVar)
                if not match:
                    continue
                wt, pos, mut = match.groups()
                pos = int(pos)
                if pos > len(protein) or protein[pos - 1] != wt:
                    continue
                mut_protein = protein[: pos - 1] + mut + protein[pos:]
                site = pos

            # In-frame indels
            elif mutType in {"In_Frame_Del", "In_Frame_Ins"}:
                parsed = parse_inframe_variant(protVar)
                if not parsed:
                    continue
                if parsed["start"] > len(protein) or parsed["end"] > len(protein):
                    continue

                if parsed["type"] == "delins":
                    mut_protein = (
                        protein[: parsed["start"] - 1]
                        + parsed["inserted"]
                        + protein[parsed["end"] :]
                    )
                    site = parsed["start"]
                elif parsed["type"] == "del":
                    mut_protein = (
                        protein[: parsed["start"] - 1] + protein[parsed["end"] :]
                    )
                    site = parsed["start"]
                elif parsed["type"] == "dup":
                    dup_seq = protein[parsed["start"] - 1 : parsed["end"]]
                    mut_protein = (
                        protein[: parsed["end"]] + dup_seq + protein[parsed["end"] :]
                    )
                    site = parsed["end"]
                elif parsed["type"] == "ins":
                    mut_protein = (
                        protein[: parsed["start"]]
                        + parsed["inserted"]
                        + protein[parsed["start"] :]
                    )
                    site = parsed["start"]
                else:
                    continue

            # Frameshift indels
            elif mutType in {"Frame_Shift_Del", "Frame_Shift_Ins"}:
                if transcriptID not in transcript2cds:
                    continue
                cds_seq = transcript2cds[transcriptID]
                parsed_hgvsc = parse_hgvsc_variant(genVar)
                if not parsed_hgvsc:
                    continue
                var_type_hgvsc, start, info = parsed_hgvsc
                mut_protein = translate_frameshift_cds(
                    cds_seq, var_type_hgvsc, start, info
                )
                if not mut_protein:
                    continue
                site = start // 3 + 1

            else:
                continue

            if not skip_digestion:
                site, peptide = digestion_func(site, mut_protein)
            else:
                peptide = mut_protein

            peptID = f"{mutType.lower()}"
            mut_pep_id = f">{peptID}|{Barcode}|{GN}|{geneID}|{transcriptID}|{protID}|{genVar}|{protVar}|peptMutPos={site}"
            output.append([mut_pep_id, peptide.split("*")[0]])

    return output


def write_fasta(
    input_file,
    output_file,
    digestion_func,
    ensp2pep,
    refseqID2prot,
    transcript2cds,
    skip_digestion,
):
    mut_pept = mutant_pept(
        input_file,
        digestion_func,
        ensp2pep,
        refseqID2prot,
        transcript2cds,
        skip_digestion,
    )
    with open(output_file, "w") as f:
        for i in mut_pept:
            f.write(i[0] + "\n")
            f.write(i[1] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Translate SNVs and Indels into mutant peptides"
    )
    parser.add_argument("-i", "--input", required=True, help="Path to input file")
    parser.add_argument("-o", "--output", required=True, help="Path to output file")
    parser.add_argument(
        "-d",
        "--digestion",
        choices=["trypsin", "lysC"],
        default="lysC",
        help="Protease used for digestion (default: lysC)",
    )
    parser.add_argument(
        "--skip_digestion",
        action="store_true",
        help="If set, skip the peptide digestion step entirely",
    )

    if "--help" in sys.argv or "-h" in sys.argv:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    ensp2pep, refseqID2prot, transcript2cds = load_proteomes()
    digestion_func = tryptic_pep if args.digestion == "trypsin" else lysC_pep
    write_fasta(
        args.input,
        args.output,
        digestion_func,
        ensp2pep,
        refseqID2prot,
        transcript2cds,
        args.skip_digestion,
    )
    print(f"Processing data from {args.input} and saving results to {args.output}")
    if args.skip_digestion:
        print("Note: digestion step was skipped.")


if __name__ == "__main__":
    main()
