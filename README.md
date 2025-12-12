# IMPACT Trial CWL Workflow

This repository contains a Common Workflow Language (CWL) implementation of the IMPACT Trial peptide discovery pipeline.

## Pipeline Overview

The workflow performs tumor-specific peptide discovery from multiple genomic alterations:
- **SNVs and Indels** from MAF files
- **Gene fusions** from STAR-Fusion and Arriba
- **Alternative splicing events** from rMATS (SE, RI, A3SS, A5SS)

### Pipeline Steps

1. **Create GTEx Matrices** - Generate junction count matrices from GTEx reference data
2. **Filter Splice Events** - Identify tumor-specific splice events relative to GTEx
3. **Translate SNVs** - Convert SNVs/Indels to peptide sequences
4. **Translate STAR-Fusion** - Generate peptides from STAR-Fusion predictions
5. **Translate Arriba Fusions** - Generate peptides from Arriba predictions
6. **Translate Splice Events** - Convert splice events to peptides
7. **Merge FASTAs** - Combine all peptide sequences
8. **Filter Against UniProt** - Remove peptides with 100% UniProt matches

## Directory Structure

```
.
├── impact-trial.cwl                      # Main workflow
├── params/
│   └── impact-trial-inputs.yml           # Example input parameters
├── tools/
│   ├── create-gtex-matrices.cwl          # Step 1 CWL tool
│   ├── filter-splice-events.cwl          # Step 2 CWL tool
│   ├── snv-translate.cwl                 # Step 3 CWL tool
│   ├── star-fusion-translate.cwl         # Step 4 CWL tool
│   ├── arriba-fusion-translate.cwl       # Step 5 CWL tool
│   ├── translate-splice-events.cwl       # Step 6 CWL tool
│   ├── merge-fasta.cwl                   # Step 7 CWL tool
│   └── filter-fasta.cwl                  # Step 8 CWL tool
├── scripts/
│   ├── create-gtex-matrices/
│   ├── splice-junction-filtering/
│   └── translation/
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
sample_id: "YOUR_SAMPLE_ID"
maf_file:
  class: File
  path: ../data/your_sample.maf
# ... (update other paths as needed)
```

### 3. Run the Workflow

Using cwltool:

```bash
cwltool impact-trial.cwl params/impact-trial-inputs.yml
```

Using toil-cwl-runner:

```bash
toil-cwl-runner impact-trial.cwl params/impact-trial-inputs.yml
```

## Outputs

The workflow produces:

- **final_filtered_fasta** - Filtered tumor-specific peptide FASTA (primary output)
- **merged_unfiltered_fasta** - All peptides before UniProt filtering
- **blast_results** - BLAST alignment results
- **gtex_matrices** - GTEx junction and target count matrices
- **tumor_specific_splice_products** - Filtered splice events TSV
- Individual peptide FASTAs from each translation step

## Input Requirements

### Required Files

1. **GTEx Reference Data**
   - rMATS merged files (SE, RI, A3SS, A5SS)
   - Histologies metadata
   - Fusion normals (STAR and Arriba)

2. **Tumor Data**
   - Somatic MAF file
   - Fusion predictions (STAR-Fusion and Arriba)
   - rMATS junction count files (SE, RI, A3SS, A5SS)
   - annoFuse filter file

3. **Reference Genomes/Proteomes**
   - GRCh38 primary assembly FASTA
   - UniProt canonical proteome

### Parameters

- `sample_id` - Sample identifier
- `min_peptide_length` - Minimum amino acid length (default: 8)
- `n_flank_bases` - Flanking nucleotides for splice junctions (default: 20)

## Docker Images

The workflow uses the following Docker images:
- `rocker/tidyverse:latest` - For R scripts
- `python:3.11` - For Python scripts
- `ncbi/blast:latest` - For BLAST filtering

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
