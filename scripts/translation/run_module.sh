#!/bin/bash

set -euo pipefail

# -------------------------------
#  run_module.sh — peptide translation pipeline
# -------------------------------

# set directory as bash script dir path
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Set up paths
root_dir=$(git rev-parse --show-toplevel)
results_dir="$SCRIPT_DIR/results"

# Default values
min_peptide_length="8"
n_flank_bases="20"

# Help page
usage() {
  cat <<EOF
Usage: bash run_module.sh [OPTIONS]

Run SNV, fusion, and splicing translation modules in sequence.

Options:
  --dna_sample_id <string>          DNA sample ID
  --rna_sample_id <string>          RNA sample ID
  --maf_file <string>               sample MAF
  --star_fusion_file <string>       sample STAR fusion file
  --arriba_fusion_file <string>     sample arriba fusion file
  --annofuse_file <string>          sample annofuse filter file
  --min_peptide_length <int>        Minimum amino acid length required to retain translated splice junction peptides (default: ${min_peptide_length})
  --n_flank_bases <int>             Number of nucleotides flanking splice junctions to include during sequence extraction (default: ${n_flank_bases})
  -h, --help                        Show this help message and exit

EOF
}

# -------------------------------
#  Parse arguments
# -------------------------------
while [ "$#" -gt 0 ]; do
  case "$1" in
    --dna_sample_id)
      [ "$#" -ge 2 ] || { echo "Error: --dna_sample_id requires a value" >&2; usage; exit 1; }
      dna_sample_id="$2"; shift 2 ;;
    --rna_sample_id)
      [ "$#" -ge 2 ] || { echo "Error: --rna_sample_id requires a value" >&2; usage; exit 1; }
      rna_sample_id="$2"; shift 2 ;;
    --maf_file)
      [ "$#" -ge 2 ] || { echo "Error: --maf_file requires a value" >&2; usage; exit 1; }
      maf_file="$2"; shift 2 ;;
    --star_fusion_file)
      [ "$#" -ge 2 ] || { echo "Error: --star_fusion_file requires a value" >&2; usage; exit 1; }
      star_fusion_file="$2"; shift 2 ;;
    --arriba_fusion_file)
      [ "$#" -ge 2 ] || { echo "Error: --arriba_fusion_file requires a value" >&2; usage; exit 1; }
      arriba_fusion_file="$2"; shift 2 ;;
    --annofuse_file)
      [ "$#" -ge 2 ] || { echo "Error: --annofuse_file requires a value" >&2; usage; exit 1; }
      annofuse_file="$2"; shift 2 ;;
    --min_peptide_length)
      [ "$#" -ge 2 ] || { echo "Error: --min_peptide_length requires a value" >&2; usage; exit 1; }
      min_peptide_length="$2"; shift 2 ;;
    --n_flank_bases)
      [ "$#" -ge 2 ] || { echo "Error: --n_flank_bases requires a value" >&2; usage; exit 1; }
      n_flank_bases="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    --) shift; break ;;               # end of options
    *)  echo "Unknown parameter: $1" >&2; usage; exit 1 ;;
  esac
done

# -------------------------------
#  Logging setup
# -------------------------------
mkdir -p "$results_dir"
LOGFILE="$results_dir/translation_$(date +%Y%m%d).log"
exec > >(tee -a "$LOGFILE") 2>&1

# -------------------------------
#  Run translation steps
# -------------------------------

echo "[$(date +%T)] Running SNV translation"
python 01-snv-translate.py \
  -i $root_dir/$maf_file \
  -o results/$dna_sample_id-snv-consensus-somatic-maf-peptides.fasta \
  --skip_digestion

echo "[$(date +%T)] Running STAR fusion translation"
python 02-star-fusion-translate.py \
    -i $root_dir/$star_fusion_file \
    -o results/$rna_sample_id-fusion-star-peptides.fasta \
    --skip_digestion \
    --gtex_normal "$root_dir/data/cfde-gtex-fusion-starfusion.tsv.gz" \
    --annofuse_filter $root_dir/$annofuse_file

echo "[$(date +%T)] Running Arriba fusion translation"
python 03-arriba-fusion-translate.py \
    -i $root_dir/$arriba_fusion_file \
    -o "results/$rna_sample_id-fusion-arriba-peptides.fasta" \
    --skip_digestion \
    --gtex_normal "$root_dir/data/cfde-gtex-arriba_formatted.annotated.tsv.gz" \
    --annofuse_filter $root_dir/$annofuse_file

# -------------------------------
#  Run alternative splicing translations
# -------------------------------

splice_event_file=$root_dir/analyses/splice-junction-filtering/results/$rna_sample_id-tumor-specific-splice-products.tsv

Rscript --vanilla 04-translate-splice-events.R --input $splice_event_file --sample_id $rna_sample_id

echo "✅ Translation completed."

# -------------------------------
#  Post-processing
# -------------------------------
bash 05-merge-fasta.sh
echo "✅ Merge completed."

bash 06-filter-fasta.sh
echo "✅ Blast completed."