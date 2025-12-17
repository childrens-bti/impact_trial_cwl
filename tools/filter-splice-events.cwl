
cwlVersion: v1.2
class: CommandLineTool

label: Filter tumor splice events for specificity
doc: |
  Filter tumor splice events junctions and targets for specificity relative to GTEx

requirements:
  DockerRequirement:
    dockerPull: "pgc-images.sbgenomics.com/childrens-bti/impact_trial_cwl:latest"
  InitialWorkDirRequirement:
    listing:
      - $(inputs.gtex_se_junction_matrix)
      - $(inputs.gtex_se_target_matrix)
      - $(inputs.gtex_ri_junction_matrix)
      - $(inputs.gtex_ri_target_matrix)
      - $(inputs.gtex_a3ss_junction_matrix)
      - $(inputs.gtex_a5ss_junction_matrix)

baseCommand: [Rscript, --vanilla, /opt/impact_trial/scripts/splice-junction-filtering/filter-splice-events.R]

arguments:
  - prefix: --tumor_se_rmats
    valueFrom: $(inputs.tumor_se_rmats.path)
  - prefix: --tumor_ri_rmats
    valueFrom: $(inputs.tumor_ri_rmats.path)
  - prefix: --tumor_a3ss_rmats
    valueFrom: $(inputs.tumor_a3ss_rmats.path)
  - prefix: --tumor_a5ss_rmats
    valueFrom: $(inputs.tumor_a5ss_rmats.path)
  - prefix: --sample_id
    valueFrom: $(inputs.sample_id)

inputs:
  sample_id:
    type: string
    doc: "Sample identifier"
  
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
  
  gtex_se_junction_matrix:
    type: File
    doc: "GTEx SE junction matrix from previous step"
  
  gtex_se_target_matrix:
    type: File
    doc: "GTEx SE target matrix from previous step"
  
  gtex_ri_junction_matrix:
    type: File
    doc: "GTEx RI junction matrix from previous step"
  
  gtex_ri_target_matrix:
    type: File
    doc: "GTEx RI target matrix from previous step"
  
  gtex_a3ss_junction_matrix:
    type: File
    doc: "GTEx A3SS junction matrix from previous step"
  
  gtex_a5ss_junction_matrix:
    type: File
    doc: "GTEx A5SS junction matrix from previous step"

outputs:
  tumor_specific_splice_products:
    type: File
    outputBinding:
      glob: "$(inputs.sample_id)-tumor-specific-splice-products.tsv"
