
cwlVersion: v1.2
class: Workflow

label: IMPACT Trial Peptide Discovery Pipeline
doc: |
  Complete workflow for tumor-specific peptide discovery from:
  - SNVs/Indels (MAF files)
  - Gene fusions (STAR-Fusion and Arriba)
  - Alternative splicing events (rMATS: SE, RI, A3SS, A5SS)
  
  Pipeline steps:
  1. Filter tumor splice events for specificity
  2. Translate SNVs to peptides
  3. Translate STAR-Fusion events to peptides
  4. Translate Arriba fusion events to peptides
  5. Translate splice events to peptides
  6. Merge all peptide FASTAs
  7. Filter against UniProt canonical proteome

requirements:
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}

inputs:
  # Sample metadata
  output_basename_dna:
    type: string
    doc: "Base name for DNA-derived output files (SNVs). Example: 'BS_1CZPVCXR_taskid'"
  output_basename_rna:
    type: string
    doc: "Base name for RNA-derived output files (fusions, splicing). Example: 'BS_S6Q7NKA3_taskid'"
  output_basename_combined:
    type: string
    doc: "Base name for combined DNA+RNA output files (merged, filtered). Example: 'BS_1CZPVCXR_BS_S6Q7NKA3_taskid'"
  
  # Pre-computed GTEx matrices
  gtex_se_junction_matrix:
    type: File
    doc: "Pre-computed GTEx SE junction matrix"
    sbg:suggestedValue:
      class: File
      path: 69441df4e0dee132ae3749ff
      name: gtex-se-junction-mat-raw.qs2
  gtex_se_target_matrix:
    type: File
    doc: "Pre-computed GTEx SE target matrix"
    sbg:suggestedValue:
      class: File
      path: 69441df3e0dee132ae3749e0
      name: gtex-se-target-mat-raw.qs2
  gtex_ri_junction_matrix:
    type: File
    doc: "Pre-computed GTEx RI junction matrix"
    sbg:suggestedValue:
      class: File
      path: 69441df3e0dee132ae3749ed
      name: gtex-ri-junction-mat-raw.qs2
  gtex_ri_target_matrix:
    type: File
    doc: "Pre-computed GTEx RI target matrix"
    sbg:suggestedValue:
      class: File
      path: 69441df4e0dee132ae374a16
      name: gtex-ri-target-mat-raw.qs2
  gtex_a3ss_junction_matrix:
    type: File
    doc: "Pre-computed GTEx A3SS junction matrix"
    sbg:suggestedValue:
      class: File
      path: 69441df4e0dee132ae374a10
      name: gtex-a3ss-junction-mat-raw.qs2
  gtex_a5ss_junction_matrix:
    type: File
    doc: "Pre-computed GTEx A5SS junction matrix"
    sbg:suggestedValue:
      class: File
      path: 69441df4e0dee132ae374a1e
      name: gtex-a5ss-junction-mat-raw.qs2
  
  gencode_gtf:
    type: File
    doc: "GENCODE gene annotation GTF file"
    sbg:suggestedValue:
      class: File
      path: 69441df4e0dee132ae374a13
      name: gencode.v39.primary_assembly.annotation.gtf.gz
  
  # Tumor splice event files
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
  
  # Proteome reference files
  ensembl_proteome:
    type: File
    doc: "Ensembl proteome reference (Homo_sapiens.GRCh38.pep.all.fa.gz)"
    sbg:suggestedValue:
      class: File
      path: 69441df4e0dee132ae3749fb
      name: Homo_sapiens.GRCh38.pep.all.fa.gz
  refseq_proteome:
    type: File
    doc: "RefSeq proteome reference (GCF_000001405.40_GRCh38.p14_protein.faa.gz)"
    sbg:suggestedValue:
      class: File
      path: 69441df4e0dee132ae374a05
      name: GCF_000001405.40_GRCh38.p14_protein.faa.gz
  cds_file:
    type: File
    doc: "CDS reference file (Homo_sapiens.GRCh38.cds.all.fa.gz)"
    sbg:suggestedValue:
      class: File
      path: 69441df3e0dee132ae3749f0
      name: Homo_sapiens.GRCh38.cds.all.fa.gz
  
  star_fusion_file:
    type: File
    doc: "STAR-Fusion predictions file"
  arriba_fusion_file:
    type: File
    doc: "Arriba fusion predictions file"
  gtex_star_normals:
    type: File
    doc: "GTEx STAR-Fusion normals"
    sbg:suggestedValue:
      class: File
      path: 69441df3e0dee132ae3749f4
      name: cfde-gtex-fusion-starfusion.tsv.gz
  gtex_arriba_normals:
    type: File
    doc: "GTEx Arriba normals"
    sbg:suggestedValue:
      class: File
      path: 69441df4e0dee132ae374a24
      name: cfde-gtex-arriba_formatted.annotated.tsv.gz
  annofuse_filter:
    type: File
    doc: "annoFuse filter file"
  genome_fasta:
    type: File
    doc: "Reference genome FASTA (GRCh38)"
    sbg:suggestedValue:
      class: File
      path: 69441df4e0dee132ae374a0b
      name: GRCh38.primary_assembly.genome.fa.gz
  
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
    sbg:suggestedValue:
      class: File
      path: 69441d765ddfaa35efbfa728
      name: UP000005640_9606.fasta

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
  
  full_matches:
    type: File
    outputSource: filter_fasta/full_matches
    doc: "List of peptide IDs with 100% UniProt matches"
  
  # Intermediate outputs
  tumor_specific_splice_products:
    type: File
    outputSource: filter_splice_events/tumor_specific_splice_products
    doc: "Filtered tumor-specific splice products"
  
  splice_junction_translation_report:
    type: File
    outputSource: translate_splice_events/junction_translation_report
    doc: "Translation report for splice junctions"
  
  splice_target_translation_report:
    type: File
    outputSource: translate_splice_events/target_translation_report
    doc: "Translation report for splice targets"
  
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
  # Step 1: Filter splice events for tumor specificity
  filter_splice_events:
    run: tools/filter-splice-events.cwl
    in:
      output_basename: output_basename_rna
      tumor_se_rmats: tumor_se_rmats
      tumor_ri_rmats: tumor_ri_rmats
      tumor_a3ss_rmats: tumor_a3ss_rmats
      tumor_a5ss_rmats: tumor_a5ss_rmats
      gtex_se_junction_matrix: gtex_se_junction_matrix
      gtex_se_target_matrix: gtex_se_target_matrix
      gtex_ri_junction_matrix: gtex_ri_junction_matrix
      gtex_ri_target_matrix: gtex_ri_target_matrix
      gtex_a3ss_junction_matrix: gtex_a3ss_junction_matrix
      gtex_a5ss_junction_matrix: gtex_a5ss_junction_matrix
      gencode_gtf: gencode_gtf
    out:
      - tumor_specific_splice_products
  
  # Step 2: Translate SNVs
  snv_translate:
    run: tools/snv-translate.cwl
    in:
      output_basename: output_basename_dna
      maf_file: maf_file
      ensembl_proteome: ensembl_proteome
      refseq_proteome: refseq_proteome
      cds_file: cds_file
    out:
      - snv_peptides
  
  # Step 3: Translate STAR-Fusion
  star_fusion_translate:
    run: tools/star-fusion-translate.cwl
    in:
      output_basename: output_basename_rna
      star_fusion_file: star_fusion_file
      gtex_normals: gtex_star_normals
      annofuse_filter: annofuse_filter
    out:
      - fusion_peptides
  
  # Step 4: Translate Arriba fusions
  arriba_fusion_translate:
    run: tools/arriba-fusion-translate.cwl
    in:
      output_basename: output_basename_rna
      arriba_fusion_file: arriba_fusion_file
      gtex_normals: gtex_arriba_normals
      annofuse_filter: annofuse_filter
    out:
      - fusion_peptides
  
  # Step 5: Translate splice events
  translate_splice_events:
    run: tools/translate-splice-events.cwl
    in:
      output_basename: output_basename_rna
      tumor_specific_splice_products: filter_splice_events/tumor_specific_splice_products
      genome_fasta: genome_fasta
      min_peptide_length: min_peptide_length
      n_flank_bases: n_flank_bases
    out:
      - splice_peptides
      - junction_translation_report
      - target_translation_report
  
  # Step 6: Merge all peptide FASTAs
  merge_fasta:
    run: tools/merge-fasta.cwl
    in:
      output_basename: output_basename_combined
      snv_peptides: snv_translate/snv_peptides
      star_fusion_peptides: star_fusion_translate/fusion_peptides
      arriba_fusion_peptides: arriba_fusion_translate/fusion_peptides
      splice_peptides: translate_splice_events/splice_peptides
    out:
      - merged_fasta
  
  # Step 7: Filter against UniProt
  filter_fasta:
    run: tools/filter-fasta.cwl
    in:
      output_basename: output_basename_combined
      custom_fasta: merge_fasta/merged_fasta
      uniprot_fasta: uniprot_fasta
    out:
      - filtered_fasta
      - blast_results
      - full_matches

$namespaces:
  sbg: "https://sevenbridges.com/"
hints:
- class: "sbg:maxNumberOfParallelInstances"
  value: 2
"sbg:links":
- id: "https://github.com/childrens-bti/impact_trial_cwl/releases/tag/v1.0.1"
  label: github-release