#!/bin/bash
set -e
set -o pipefail

if [ -f .env ]; then
  source .env
fi

PROFILE_NAME="impact-sso"

echo "🔍 Checking AWS SSO session for profile: $PROFILE_NAME"

if ! aws sts get-caller-identity --profile "$PROFILE_NAME" > /dev/null 2>&1; then
  echo "🔐 AWS SSO session expired or missing. Logging in..."
  aws sso login --profile "$PROFILE_NAME" || {
    echo "❌ AWS SSO login failed. Please check your credentials."
    exit 1
  }
else
  echo "✅ AWS SSO session is active."
fi

# Define S3 URI and local release path
RELEASE="v8"
URI="s3://bti-private-us-east-1-prd-impact-trial/testing"
DEST="data/$RELEASE"

# Remove old symlinks in data
find data -type l -delete

mkdir -p "$DEST"

MD5_FILE="$DEST/md5sum.txt"
aws s3 cp "$URI/$RELEASE/md5sum.txt" "$MD5_FILE" --profile "$PROFILE_NAME"

# Consider the filenames in the md5sum file and the release notes
FILES=(`tr -s ' ' < data/$RELEASE/md5sum.txt | cut -d ' ' -f 2` release-notes.md)

echo "📦 Downloading files listed in md5sum.txt..."
for file in "${FILES[@]}"; do
  # skip empty lines defensively
  [[ -z "$file" ]] && continue
  src="$URI/$RELEASE/$file"
  dst="$DEST/$file"
  aws s3 cp "$src" "$dst" --profile "$PROFILE_NAME"
done

# External GENCODE/Ensembl/RefSeq files
cd "data"
# 1) GENCODE v39 annotations 
GENCODE_ANNOT="https://bti-openaccess-us-east-1-prd-references.s3.us-east-1.amazonaws.com/gencode.v39.primary_assembly.annotation.gtf.gz"
# 2) Ensembl proteome
ENS_PROT="https://bti-openaccess-us-east-1-prd-references.s3.us-east-1.amazonaws.com/Homo_sapiens.GRCh38.pep.all.fa.gz"
# 3) Refseq proteome
REFSEQ_PROT="https://bti-openaccess-us-east-1-prd-references.s3.us-east-1.amazonaws.com/GCF_000001405.40_GRCh38.p14_protein.faa.gz"
# 4) GENCODE v39 FASTA
GENCODE_FASTA="https://bti-openaccess-us-east-1-prd-references.s3.us-east-1.amazonaws.com/GRCh38.primary_assembly.genome.fa.gz"
# 5) Ensembl cds
ENS_CDS="https://bti-openaccess-us-east-1-prd-references.s3.us-east-1.amazonaws.com/Homo_sapiens.GRCh38.cds.all.fa.gz"

if [ ! -e "${GENCODE_ANNOT##*/}" ]; then
  echo "Downloading ${GENCODE_ANNOT##*/}"
  wget -nc "$GENCODE_ANNOT"
fi

if [ ! -e "${ENS_PROT##*/}" ]; then
  echo "Downloading ${ENS_PROT##*/}"
  wget -nc "$ENS_PROT"
fi

if [ ! -e "${REFSEQ_PROT##*/}" ]; then
  echo "Downloading ${REFSEQ_PROT##*/}"
  wget -nc "$REFSEQ_PROT"
fi

if [ ! -e "${GENCODE_FASTA##*/}" ]; then
  echo "Downloading ${GENCODE_FASTA##*/}"
  wget -nc "$GENCODE_FASTA"
fi

if [ ! -e "${ENS_CDS##*/}" ]; then
  echo "Downloading ${ENS_CDS##*/}"
  wget -nc "$ENS_CDS"
fi

# GTEx fusion files
gtex_arriba="https://bti-openaccess-us-east-1-prd-rokita-lab.s3.us-east-1.amazonaws.com/data-merges/v1/cfde-gtex-arriba_formatted.annotated.tsv.gz"
gtex_star="https://bti-openaccess-us-east-1-prd-rokita-lab.s3.us-east-1.amazonaws.com/data-merges/v1/cfde-gtex-fusion-starfusion.tsv.gz"

if [ ! -e "${gtex_arriba##*/}" ]; then
  echo "Downloading ${gtex_arriba##*/}"
  wget -nc "$gtex_arriba"
fi

if [ ! -e "${gtex_star##*/}" ]; then
  echo "Downloading ${gtex_star##*/}"
  wget -nc "$gtex_star"
fi


# Validate MD5s
cd "$RELEASE"
echo "Checking MD5 hashes..."
md5sum -c md5sum.txt || echo "⚠️ Warning: some checksums failed"
cd ../..

# Symlink files from release to top-level data dir
for file in "${FILES[@]}"; do
  # Skip if file is empty or already exists as a directory
  [[ -z "$file" || "$file" == "$RELEASE" ]] && continue
  ln -sfn "$RELEASE/$file" "data/$file"
done

echo "✅ All done. Files available under data/ and data/$RELEASE/"