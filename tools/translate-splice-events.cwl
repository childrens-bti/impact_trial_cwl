
cwlVersion: v1.2
class: CommandLineTool

label: Translate splice events to peptides
doc: |
  Performs in silico translation of tumor-specific splice products

requirements:
  DockerRequirement:
    dockerPull: "pgc-images.sbgenomics.com/childrens-bti/impact_trial_cwl:latest"
  InitialWorkDirRequirement:
    listing:
      - $(inputs.genome_fasta)

baseCommand: [Rscript, --vanilla, /opt/impact_trial/scripts/translation/04-translate-splice-events.R]

arguments:
  - prefix: --input
    valueFrom: $(inputs.tumor_specific_splice_products.path)
  - prefix: --sample_id
    valueFrom: $(inputs.output_basename)
  - prefix: --min_peptide_length
    valueFrom: $(inputs.min_peptide_length)
  - prefix: --n_flank_bases
    valueFrom: $(inputs.n_flank_bases)

inputs:
  output_basename:
    type: string
    doc: "Base name for output files (passed as sample_id to script)"
  
  tumor_specific_splice_products:
    type: File
    doc: "Tumor-specific splice products TSV from filter step"
  
  genome_fasta:
    type: File
    doc: "Reference genome FASTA"
  
  min_peptide_length:
    type: int
    default: 8
    doc: "Minimum peptide length"
  
  n_flank_bases:
    type: int
    default: 20
    doc: "Number of junction-flanking nucleotides"

outputs:
  splice_peptides:
    type: File
    outputBinding:
      glob: "$(inputs.output_basename)_tumor_specific_splice_event_peptides.fasta"
  
  junction_translation_report:
    type: File
    outputBinding:
      glob: "$(inputs.output_basename)-tumor-specific-splice-junction-translation-report.tsv"
  
  target_translation_report:
    type: File
    outputBinding:
      glob: "$(inputs.output_basename)-tumor-specific-splice-target-translation-report.tsv"
