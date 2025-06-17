# SNV and Fusion Peptide Translation Module

## Usage

`bash run_module.sh`

## Folder contents

1. `01-snv-translate.py`: creates a fasta file consists of peptide sequences with SNVs called from WES.
2. `02-start-fusion-translate.py`: create a fasta file that has fusion protein sequences from the STAR fusion output.
3. `03-arriba-fusion-translate.py`: create a fasta file that has fusion protein sequences from the Arriba fusion output.

### General usage of scripts

#### `01-snv-translate.py`

**script description**
```
python 01-snv-translate.py -h

usage: 01-snv-translate.py [-h] -i INPUT -o OUTPUT -d {trypsin,lysC}

A python script for SNVs

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to input file
  -o OUTPUT, --output OUTPUT
                        Path to output file
  -d {trypsin,lysC}, --digestion {trypsin,lysC}
                        Enzymes: trypsin or lysC
```

```
# Run script
python 01-snv-translate.py -i ../data/snv-consensus-plus-hotspots.maf-test.tsv -o snv-consensus-plus-hotspots.maf.fasta -d lysC
```

#### `02-star-fusion-translate.py`

**script description**
```
python 02-star-fusion-translate.py -h

usage: 02-star-fusion-translate.py [-h] -i INPUT -o OUTPUT [-d {trypsin,lysC}]

A python script for star fusions

options:
  -h, --help           show this help message and exit
  -i, --input INPUT    Path to fusion TSV file
  -o, --output OUTPUT  Path to output FASTA file
  -d, --digestion {trypsin,lysC}
                        Protease used for digestion (default: lysC)
```

```
# Run script
python 02-star-fusion-translate.py -i ../data/fusion-starfusion-test.tsv.gz -o fusion-starfusion-peptides.fasta
```
Fusion sequences are validated by concatenating the left and right CDS regions and translating them into protein sequences using the Homo_sapiens.GRCh38.cds.all.fa reference downloaded from ensembl ftp site: https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz

#### `03-arriba-fusion-translate.py`

**script description**
```
python 03-arriba-fusion-translate.py -h

usage: 03-arriba-fusion-translate.py [-h] -i INPUT -o OUTPUT [-d {trypsin,lysC}]

A python script for star fusions

options:
  -h, --help           show this help message and exit
  -i, --input INPUT    Path to fusion TSV file
  -o, --output OUTPUT  Path to output FASTA file
  -d, --digestion {trypsin,lysC}
                        Protease used for digestion (default: lysC)
```

```
# Run script
python 03-arriba-fusion-translate.py -i ../data/fusion-starfusion-test.tsv.gz -o fusion-starfusion-peptides.fasta
```

## Analysis module directory structure

```
.
├── 01-snv-translate.py
├── 02-star-fusion-translate.py
├── 03-arriba-fusion-translate.py
├── 02-pca-mirna-batch-comparison.html
├── README.md
├── results
│   ├── snv-consensus-plus-hotspots.maf.fasta
│   ├── fusion-starfusion-peptides.fasta
│   └── fusion-arriba-peptides.fasta
└── run_module.sh
```