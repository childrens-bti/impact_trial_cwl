# Unified Peptide Translation Module

Generate peptide FASTA files from:
- SNVs and Indels
- Fusions (STAR-Fusion and/or Arriba)
- Alternative splicing events (rMATS: SE, RI, A3SS, A5SS)

## Usage

`bash run_module.sh`

## Folder contents

1. `01-snv-translate.py`: creates a fasta file consists of peptide sequences with SNVs and Indels called from WES.
2. `02-star-fusion-translate.py`: create a fasta file that has fusion protein sequences from the STAR fusion output.
3. `03-arriba-fusion-translate.py`: create a fasta file that has fusion protein sequences from the Arriba fusion output.
4. `04-rMATS_SE_translation.R`: translates skipped exon (SE) splice events
5. `05-rMATS_RI_translation.R`: translates retained intron (RI) splice events
6. `06-rMATS_A3SS_translation.R`: translates retained alternative 3' splice site (A3SS) usage events
7. `07-rMATS_A5SS_translation.R`: translates retained alternative 5' splice site (A5SS) usage events
8. `08-merge-fasta.sh`: merges all results into one custom fasta file
9. `09-filter-fasta.sh`: blasts `custom.fasta` peptides against the uniprot canonical peptide BLAST database (UP000005640_9606.fasta, obtained from UniProt, last modified date of 2024-07-24), and filters out any peptides with 100% query match. 

### General usage of scripts

#### `01-snv-translate.py`

**script description**
```
python 01-snv-translate.py -h

usage: 01-snv-translate.py [-h] -i INPUT -o OUTPUT [-d {trypsin,lysC}] [--skip_digestion]

Translate SNVs and Indels into mutant peptides

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     Path to input file
  -o, --output OUTPUT   Path to output file
  -d {trypsin,lysC}, --digestion {trypsin,lysC}
                        Protease used for digestion (default: lysC)
  --skip_digestion      If set, skip the peptide digestion step entirely
```

#### `02-star-fusion-translate.py`

**script description**
```
python 02-star-fusion-translate.py -h

usage: 02-star-fusion-translate.py [-h] -i INPUT -o OUTPUT [-d {trypsin,lysC}] [--skip_digestion] [--gtex_normals GTEX_NORMALS] [--annofuse_filter ANNOFUSE_FILTER]

A python script for STAR fusions

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to fusion TSV file
  -o OUTPUT, --output OUTPUT
                        Path to output FASTA file
  -d {trypsin,lysC}, --digestion {trypsin,lysC}
                        Protease used for digestion (default: lysC)
  --skip_digestion      If set, skip the peptide digestion step entirely
  --gtex_normals GTEX_NORMALS
                        Path to GTEx normals TSV with column 'FusionName'
  --annofuse_filter ANNOFUSE_FILTER
                        Path to *.annoFuse_filter.tsv (must contain 'FusionName'); only these fusions will be translated
```

Fusion sequences are validated by concatenating the left and right CDS regions and translating them into protein sequences using the Homo_sapiens.GRCh38.cds.all.fa reference downloaded from ensembl ftp site: https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz

#### `03-arriba-fusion-translate.py`

**script description**
```
python 03-arriba-fusion-translate.py -h

usage: 03-arriba-fusion-translate.py [-h] -i INPUT -o OUTPUT [-d {trypsin,lysC}] [--skip_digestion] [--gtex_normals GTEX_NORMALS] [--annofuse_filter ANNOFUSE_FILTER]

A python script for Arriba fusions

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to fusion TSV file
  -o OUTPUT, --output OUTPUT
                        Path to output FASTA file
  -d {trypsin,lysC}, --digestion {trypsin,lysC}
                        Protease used for digestion (default: lysC)
  --skip_digestion      If set, skip the peptide digestion step entirely
  --gtex_normals GTEX_NORMALS
                        Path to GTEx normals TSV with an ordered column 'gene1--gene2'
  --annofuse_filter ANNOFUSE_FILTER
                        Path to *.annoFuse_filter.tsv (must contain 'FusionName'); only these fusions will be translated
```

## Analysis module directory structure

```
.
├── 01-snv-translate.py
├── 02-star-fusion-translate.py
├── 03-arriba-fusion-translate.py
├── 04-rMATS_SE_translation.R
├── 05-rMATS_RI_translation.R
├── 06-rMATS_A3SS_translation.R
├── 07-rMATS_A5SS_translation.R
├── 08-merge-fasta.sh
├── 09-filter-fasta.sh
├── input
│   └── UP000005640_9606.fasta
├── README.md
├── results
│   ├── A3SS-event-translation-report.tsv
│   ├── A5SS-event-translation-report.tsv
│   ├── blast_results.tsv
│   ├── custom_filtered.fasta
│   ├── custom.fasta
│   ├── full_matches.txt
│   ├── fusion-arriba-peptides.fasta
│   ├── RI-event-translation-report.tsv
│   ├── rMATS_A3SS_peptides.fasta
│   ├── rMATS_A5SS_peptides.fasta
│   ├── rMATS_RI_peptides.fasta
│   ├── rMATS_SE_peptides.fasta
│   ├── SE-event-translation-report.tsv
│   ├── snv-consensus-somatic-maf-peptides.fasta
│   ├── translation_20250826.log
│   ├── uniprot_db.pdb
│   ├── uniprot_db.phr
│   ├── uniprot_db.pin
│   ├── uniprot_db.pjs
│   ├── uniprot_db.pot
│   ├── uniprot_db.psq
│   ├── uniprot_db.ptf
│   └── uniprot_db.pto
├── run_module.sh
└── util
    ├── annofuse_filtering.py
    └── digestion.py
```