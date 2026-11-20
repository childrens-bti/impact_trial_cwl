
cwlVersion: v1.2
class: CommandLineTool

label: Translate SNVs and Indels to peptides
doc: |
  Creates a FASTA file consisting of peptide sequences with SNVs and Indels called from WES

requirements:
  DockerRequirement:
    dockerPull: "pgc-images.sbgenomics.com/childrens-bti/impact_trial_cwl:v1.3.0"

baseCommand: [python, /opt/impact_trial/scripts/translation/01-snv-translate.py]

arguments:
  - prefix: -i
    valueFrom: $(inputs.maf_file.path)
  - prefix: -o
    valueFrom: $(inputs.output_basename)_snv-peptides.fasta
  - prefix: --ensembl_proteome
    valueFrom: $(inputs.ensembl_proteome.path)
  - prefix: --refseq_proteome
    valueFrom: $(inputs.refseq_proteome.path)
  - prefix: --cds_file
    valueFrom: $(inputs.cds_file.path)
  - --skip_digestion

inputs:
  output_basename:
    type: string
    doc: "Base name for output files"
  
  maf_file:
    type: File
    doc: "Somatic MAF file"
  
  ensembl_proteome:
    type: File
    doc: "Ensembl proteome reference file"
  
  refseq_proteome:
    type: File
    doc: "RefSeq proteome reference file"
  
  cds_file:
    type: File
    doc: "CDS reference file"

outputs:
  snv_peptides:
    type: File
    outputBinding:
      glob: "$(inputs.output_basename)_snv-peptides.fasta"
