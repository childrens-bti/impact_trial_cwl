
cwlVersion: v1.2
class: CommandLineTool

label: Translate STAR-Fusion to peptides
doc: |
  Creates a FASTA file with fusion protein sequences from STAR-Fusion output

requirements:
  DockerRequirement:
    dockerPull: "pgc-images.sbgenomics.com/childrens-bti/impact_trial_cwl:latest"

baseCommand: [python, /opt/impact_trial/scripts/translation/02-star-fusion-translate.py]

arguments:
  - prefix: -i
    valueFrom: $(inputs.star_fusion_file.path)
  - prefix: -o
    valueFrom: $(inputs.sample_id)_star-fusion-peptides.fasta
  - --skip_digestion
  - prefix: --gtex_normals
    valueFrom: $(inputs.gtex_normals.path)
  - prefix: --annofuse_filter
    valueFrom: $(inputs.annofuse_filter.path)

inputs:
  sample_id:
    type: string
    doc: "Sample identifier for output filename"
  
  star_fusion_file:
    type: File
    doc: "STAR-Fusion predictions file"
  
  gtex_normals:
    type: File
    doc: "GTEx normal fusion file"
  
  annofuse_filter:
    type: File
    doc: "annoFuse filter file"

outputs:
  fusion_peptides:
    type: File
    outputBinding:
      glob: "$(inputs.sample_id)_star-fusion-peptides.fasta"
