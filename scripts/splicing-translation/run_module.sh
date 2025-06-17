#!/bin/bash

set -e
set -o pipefail

# Run DE analysis script
Rscript -e "rmarkdown::render('rMATS_SE_translation.Rmd')"