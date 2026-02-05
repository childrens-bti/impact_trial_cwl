cwlVersion: v1.2
class: CommandLineTool

label: Merge FASTA files
doc: |
  Concatenates all peptide FASTA files into a single merged FASTA

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.snv_peptides)
      - $(inputs.star_fusion_peptides)
      - $(inputs.arriba_fusion_peptides)
      - $(inputs.splice_peptides)

baseCommand: [cat]

inputs:
  output_basename:
    type: string
    doc: "Base name for output files"
  
  snv_peptides:
    type: File
    inputBinding:
      position: 1
    doc: "SNV peptides FASTA"
  
  star_fusion_peptides:
    type: File
    inputBinding:
      position: 2
    doc: "STAR-Fusion peptides FASTA"
  
  arriba_fusion_peptides:
    type: File
    inputBinding:
      position: 3
    doc: "Arriba fusion peptides FASTA"
  
  splice_peptides:
    type: File
    inputBinding:
      position: 4
    doc: "Splice event peptides FASTA"

stdout: $(inputs.output_basename)_all_peptides_merged.fasta

outputs:
  merged_fasta:
    type: stdout
