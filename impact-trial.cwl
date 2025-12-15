#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow

label: IMPACT Trial Peptide Discovery Pipeline
doc: |
  Complete workflow for tumor-specific peptide discovery from:
  - SNVs/Indels (MAF files)
  - Gene fusions (STAR-Fusion and Arriba)
  - Alternative splicing events (rMATS: SE, RI, A3SS, A5SS)
  
  Pipeline steps:
  1. Create GTEx junction count matrices
  2. Filter tumor splice events for specificity
  3. Translate SNVs to peptides
  4. Translate STAR-Fusion events to peptides
  5. Translate Arriba fusion events to peptides
  6. Translate splice events to peptides
  7. Merge all peptide FASTAs
  8. Filter against UniProt canonical proteome

requirements:
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  # Sample metadata
  sample_id:
    type: string
    doc: "Sample identifier"
  
  # GTEx reference files (Step 1)
  gtex_se_rmats:
    type: File
    doc: "GTEx rMATS SE merged file"
  gtex_ri_rmats:
    type: File
    doc: "GTEx rMATS RI merged file"
  gtex_a3ss_rmats:
    type: File
    doc: "GTEx rMATS A3SS merged file"
  gtex_a5ss_rmats:
    type: File
    doc: "GTEx rMATS A5SS merged file"
  histologies_file:
    type: File
    doc: "GTEx histologies metadata"
  
  # Tumor splice event files (Step 2)
  tumor_se_rmats:
    type: File
    doc: "Tumor rMATS SE JC file"
  tumor_ri_rmats:
    type: File
    doc: "Tumor rMATS RI JC file"
  tumor_a3ss_rmats:
    type: File
    doc: "Tumor rMATS A3SS JC file"
  tumor_a5ss_rmats:
    type: File
    doc: "Tumor rMATS A5SS JC file"
  
  # Translation inputs (Steps 3-6)
  maf_file:
    type: File
    doc: "Somatic MAF file for SNV translation"
  star_fusion_file:
    type: File
    doc: "STAR-Fusion predictions file"
  arriba_fusion_file:
    type: File
    doc: "Arriba fusion predictions file"
  gtex_star_normals:
    type: File
    doc: "GTEx STAR-Fusion normals"
  gtex_arriba_normals:
    type: File
    doc: "GTEx Arriba normals"
  annofuse_filter:
    type: File
    doc: "annoFuse filter file"
  genome_fasta:
    type: File
    doc: "Reference genome FASTA (GRCh38)"
  
  # Translation parameters
  min_peptide_length:
    type: int
    default: 8
    doc: "Minimum peptide length for splice events"
  n_flank_bases:
    type: int
    default: 20
    doc: "Number of junction-flanking nucleotides"
  
  # Filtering input (Step 8)
  uniprot_fasta:
    type: File
    doc: "UniProt canonical proteome FASTA"
  
  # Script files
  create_gtex_matrices_script:
    type: File
    doc: "R script for creating GTEx matrices"
  filter_splice_events_script:
    type: File
    doc: "R script for filtering splice events"
  snv_translate_script:
    type: File
    doc: "Python script for SNV translation"
  star_fusion_translate_script:
    type: File
    doc: "Python script for STAR-Fusion translation"
  arriba_fusion_translate_script:
    type: File
    doc: "Python script for Arriba fusion translation"
  translate_splice_events_script:
    type: File
    doc: "R script for splice event translation"

outputs:
  # Final outputs
  final_filtered_fasta:
    type: File
    outputSource: filter_fasta/filtered_fasta
    doc: "Final filtered tumor-specific peptide FASTA"
  
  merged_unfiltered_fasta:
    type: File
    outputSource: merge_fasta/merged_fasta
    doc: "Merged unfiltered peptide FASTA"
  
  blast_results:
    type: File
    outputSource: filter_fasta/blast_results
    doc: "BLAST results against UniProt"
  
  # Intermediate outputs
  gtex_matrices:
    type:
      type: array
      items: File
    outputSource:
      - create_gtex_matrices/gtex_se_junction_matrix
      - create_gtex_matrices/gtex_se_target_matrix
      - create_gtex_matrices/gtex_ri_junction_matrix
      - create_gtex_matrices/gtex_ri_target_matrix
      - create_gtex_matrices/gtex_a3ss_junction_matrix
      - create_gtex_matrices/gtex_a5ss_junction_matrix
    doc: "GTEx junction and target matrices"
  
  tumor_specific_splice_products:
    type: File
    outputSource: filter_splice_events/tumor_specific_splice_products
    doc: "Filtered tumor-specific splice products"
  
  snv_peptides:
    type: File
    outputSource: snv_translate/snv_peptides
    doc: "SNV-derived peptides"
  
  star_fusion_peptides:
    type: File
    outputSource: star_fusion_translate/fusion_peptides
    doc: "STAR-Fusion-derived peptides"
  
  arriba_fusion_peptides:
    type: File
    outputSource: arriba_fusion_translate/fusion_peptides
    doc: "Arriba fusion-derived peptides"
  
  splice_peptides:
    type: File
    outputSource: translate_splice_events/splice_peptides
    doc: "Splice event-derived peptides"

steps:
  # Step 1: Create GTEx matrices
  create_gtex_matrices:
    run: tools/create-gtex-matrices.cwl
    in:
      script: create_gtex_matrices_script
      gtex_se_rmats: gtex_se_rmats
      gtex_ri_rmats: gtex_ri_rmats
      gtex_a3ss_rmats: gtex_a3ss_rmats
      gtex_a5ss_rmats: gtex_a5ss_rmats
      histologies_file: histologies_file
    out:
      - gtex_se_junction_matrix
      - gtex_se_target_matrix
      - gtex_ri_junction_matrix
      - gtex_ri_target_matrix
      - gtex_a3ss_junction_matrix
      - gtex_a5ss_junction_matrix
  
  # Step 2: Filter splice events for tumor specificity
  filter_splice_events:
    run: tools/filter-splice-events.cwl
    in:
      script: filter_splice_events_script
      sample_id: sample_id
      tumor_se_rmats: tumor_se_rmats
      tumor_ri_rmats: tumor_ri_rmats
      tumor_a3ss_rmats: tumor_a3ss_rmats
      tumor_a5ss_rmats: tumor_a5ss_rmats
      gtex_se_junction_matrix: create_gtex_matrices/gtex_se_junction_matrix
      gtex_se_target_matrix: create_gtex_matrices/gtex_se_target_matrix
      gtex_ri_junction_matrix: create_gtex_matrices/gtex_ri_junction_matrix
      gtex_ri_target_matrix: create_gtex_matrices/gtex_ri_target_matrix
      gtex_a3ss_junction_matrix: create_gtex_matrices/gtex_a3ss_junction_matrix
      gtex_a5ss_junction_matrix: create_gtex_matrices/gtex_a5ss_junction_matrix
    out:
      - tumor_specific_splice_products
  
  # Step 3: Translate SNVs
  snv_translate:
    run: tools/snv-translate.cwl
    in:
      script: snv_translate_script
      maf_file: maf_file
    out:
      - snv_peptides
  
  # Step 4: Translate STAR-Fusion
  star_fusion_translate:
    run: tools/star-fusion-translate.cwl
    in:
      script: star_fusion_translate_script
      star_fusion_file: star_fusion_file
      gtex_normals: gtex_star_normals
      annofuse_filter: annofuse_filter
    out:
      - fusion_peptides
  
  # Step 5: Translate Arriba fusions
  arriba_fusion_translate:
    run: tools/arriba-fusion-translate.cwl
    in:
      script: arriba_fusion_translate_script
      arriba_fusion_file: arriba_fusion_file
      gtex_normals: gtex_arriba_normals
      annofuse_filter: annofuse_filter
    out:
      - fusion_peptides
  
  # Step 6: Translate splice events
  translate_splice_events:
    run: tools/translate-splice-events.cwl
    in:
      script: translate_splice_events_script
      sample_id: sample_id
      tumor_specific_splice_products: filter_splice_events/tumor_specific_splice_products
      genome_fasta: genome_fasta
      min_peptide_length: min_peptide_length
      n_flank_bases: n_flank_bases
    out:
      - splice_peptides
  
  # Step 7: Merge all peptide FASTAs
  merge_fasta:
    run: tools/merge-fasta.cwl
    in:
      snv_peptides: snv_translate/snv_peptides
      star_fusion_peptides: star_fusion_translate/fusion_peptides
      arriba_fusion_peptides: arriba_fusion_translate/fusion_peptides
      splice_peptides: translate_splice_events/splice_peptides
    out:
      - merged_fasta
  
  # Step 8: Filter against UniProt
  filter_fasta:
    run: tools/filter-fasta.cwl
    in:
      custom_fasta: merge_fasta/merged_fasta
      uniprot_fasta: uniprot_fasta
    out:
      - filtered_fasta
      - blast_results
