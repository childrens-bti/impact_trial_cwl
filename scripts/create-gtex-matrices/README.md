# Unified Peptide Translation Module

Generate GTEx splice event junction and target matrices:

## Usage

`bash run_module.sh`

## Folder contents

1. `create-gtex-jc-count-matrices.R`: Converts GTEx rMATS splice event output to mean junction count matrices by GTEx subgroup

### General usage of scripts

#### `bash run_module.sh`

## Analysis module directory structure

```
.
├── README.md
├── create-gtex-jc-count-matrices.R
└── results
    ├── gtex-a3ss-junction-mat-raw.qs2
    ├── gtex-a5ss-junction-mat-raw.qs2
    ├── gtex-ri-junction-mat-raw.qs2
    ├── gtex-ri-target-mat-raw.qs2
    ├── gtex-se-junction-mat-raw.qs2
    └── gtex-se-target-mat-raw.qs2
```