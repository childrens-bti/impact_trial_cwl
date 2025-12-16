# Splice event filtering module

Filter tumor splice events junctions and targets for specificity relative to GTEx 

## Usage

```
Rscript --vanilla filter-splice-events.R \
--tumor_se_rmats <rmats_se_file> \
--tumor_ri_rmats <rmats_ri_file> \
--tumor_a3ss_rmats <rmats_a3ss_file> \
--tumor_a5ss_rmats <rmats_a5ss_file>
```

## Folder contents

1. `filter-splice-events.R`: Filter junctions and targets in tumor rMATS output for tumor-specific splice products

## Analysis module directory structure

```
.
├── README.md
├── filter-splice-events.R
└── results
    └── tumor-specific-splice-products.tsv
```