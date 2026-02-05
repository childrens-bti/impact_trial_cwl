# IMPACT Trial CWL Workflow

**Latest Release:** v1.0.0 (Stable)

This repository contains a Common Workflow Language (CWL) implementation of the IMPACT Trial peptide discovery pipeline.

## Pipeline Overview

The workflow performs tumor-specific peptide discovery from multiple genomic alterations:
- **SNVs and Indels** from MAF files
- **Gene fusions** from STAR-Fusion and Arriba
- **Alternative splicing events** from rMATS (SE, RI, A3SS, A5SS)

### Pipeline Steps

1. **Filter Splice Events** - Identify tumor-specific splice events relative to pre-computed GTEx matrices
2. **Translate SNVs** - Convert SNVs/Indels to peptide sequences
3. **Translate STAR-Fusion** - Generate peptides from STAR-Fusion predictions
4. **Translate Arriba Fusions** - Generate peptides from Arriba predictions
5. **Translate Splice Events** - Convert splice events to peptides
6. **Merge FASTAs** - Combine all peptide sequences
7. **Filter Against UniProt** - Remove peptides with 100% UniProt matches

**Note:** This workflow uses pre-computed GTEx junction and target matrices to avoid the memory-intensive matrix generation step (previously Step 1, requiring 65-70GB RAM).

## Directory Structure

```
.
├── impact-trial.cwl                      # Main workflow (7 steps)
├── params/
│   └── impact-trial-inputs.yml           # Example input parameters
├── tools/
│   ├── filter-splice-events.cwl          # Step 1 CWL tool
│   ├── snv-translate.cwl                 # Step 2 CWL tool
│   ├── star-fusion-translate.cwl         # Step 3 CWL tool
│   ├── arriba-fusion-translate.cwl       # Step 4 CWL tool
│   ├── translate-splice-events.cwl       # Step 5 CWL tool
│   ├── merge-fasta.cwl                   # Step 6 CWL tool
│   └── filter-fasta.cwl                  # Step 7 CWL tool
├── scripts/
│   ├── splice-junction-filtering/
│   │   └── filter-splice-events.R
│   └── translation/
│       ├── 01-snv-translate.py
│       ├── 02-star-fusion-translate.py
│       ├── 03-arriba-fusion-translate.py
│       ├── 04-translate-splice-events.R
│       ├── 05-merge-fasta.sh
│       └── 06-filter-fasta.sh
├── data/                                  # Input data files
└── outputs/                               # Workflow outputs
```

## Requirements

- CWL runner (e.g., cwltool, toil-cwl-runner, or cromwell)
- Docker (for containerized execution)
- Input data files (see `data/` directory)

## Usage

### 1. Download Required Data

First, download all required reference and input data:

```bash
bash download_data.sh
```

### 2. Prepare Input Parameters

Edit `params/impact-trial-inputs.yml` to specify your input files and parameters:

```yaml
# Sample identifiers
dna_sample_id: "YOUR_DNA_SAMPLE_ID"  # For MAF/SNV outputs
rna_sample_id: "YOUR_RNA_SAMPLE_ID"  # For fusion/splice outputs

# Pre-computed GTEx matrices (required)
se_gtex_junction_mat:
  class: File
  path: ../data/gtex-se-junction-mat-raw.qs2
se_gtex_target_mat:
  class: File
  path: ../data/gtex-se-target-mat-raw.qs2
# ... (see params/impact-trial-inputs.yml for full example)

# Tumor data files
maf_file:
  class: File
  path: ../data/your_sample.maf
# ... (update other paths as needed)
```

### 3. Run the Workflow

Using cwltool:

```bash
cwltool --outdir outputs/patient_1 impact-trial.cwl params/impact-trial-inputs.yml
```

With additional options for temporary file management:

```bash
cwltool --leave-tmpdir --tmpdir-prefix ./.cwl-tmp/ --tmp-outdir-prefix ./.cwl-out/ \
    --outdir outputs/patient_1 impact-trial.cwl params/impact-trial-inputs.yml
```

## Outputs

The workflow produces:

### Primary Outputs
- **final_filtered_fasta** - Filtered tumor-specific peptide FASTA (peptides unique vs UniProt)
- **merged_unfiltered_fasta** - All peptides before UniProt filtering
- **blast_results** - BLAST alignment results against UniProt
- **full_matches** - List of peptide IDs with 100% UniProt matches (filtered out)

### Intermediate Outputs
- **tumor_specific_splice_products** - Filtered tumor-specific splice events TSV
- **splice_junction_translation_report** - Translation report for splice junctions
- **splice_target_translation_report** - Translation report for splice targets

All outputs are written directly to the specified output directory (e.g., `outputs/patient_1/`).

## Input Requirements

### Required Files

1. **Pre-computed GTEx Matrices** (6 files)
   - SE junction matrix: `gtex-se-junction-mat-raw.qs2`
   - SE target matrix: `gtex-se-target-mat-raw.qs2`
   - RI junction matrix: `gtex-ri-junction-mat-raw.qs2`
   - RI target matrix: `gtex-ri-target-mat-raw.qs2`
   - A3SS junction matrix: `gtex-a3ss-junction-mat-raw.qs2`
   - A5SS junction matrix: `gtex-a5ss-junction-mat-raw.qs2`

2. **GTEx Reference Data**
   - Fusion normals (STAR and Arriba)
   - annoFuse filter file

3. **Tumor Data**
   - Somatic MAF file
   - Fusion predictions (STAR-Fusion and Arriba)
   - rMATS junction count files (SE, RI, A3SS, A5SS)

4. **Reference Genomes/Proteomes**
   - GRCh38 primary assembly FASTA
   - GENCODE GTF annotation (e.g., gencode.v39.primary_assembly.annotation.gtf.gz)
   - Ensembl proteome FASTA
   - RefSeq proteome FASTA
   - CDS reference FASTA
   - UniProt canonical proteome FASTA

### Parameters

- `dna_sample_id` - DNA sample identifier (for MAF/SNV outputs)
- `rna_sample_id` - RNA sample identifier (for fusion/splice outputs)
- `min_peptide_length` - Minimum amino acid length (default: 8)
- `n_flank_bases` - Flanking nucleotides for splice junctions (default: 20)

## Docker Image

The workflow uses a unified Docker image containing all required dependencies:
- **Image:** `pgc-images.sbgenomics.com/childrens-bti/impact_trial_cwl:latest`
- **Base:** `rocker/tidyverse:4.4.0`
- **Includes:** R packages (tidyverse, qs2, Biostrings, rtracklayer, optparse), Python 3.11, rpg 2.0.3, NCBI BLAST+, bedtools, samtools, seqkit

## Development

### Adding New Tools

1. Create a new CWL CommandLineTool in `tools/`
2. Update the main workflow in `impact-trial.cwl`
3. Add corresponding script files in `scripts/`

### Testing

Test individual tools:

```bash
cwltool tools/snv-translate.cwl test-inputs.yml
```

Validate workflow:

```bash
cwltool --validate impact-trial.cwl
```

## License

See LICENSE file for details.
