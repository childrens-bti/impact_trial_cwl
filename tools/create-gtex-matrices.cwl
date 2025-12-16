
cwlVersion: v1.2
class: CommandLineTool

label: Create GTEx junction count matrices
doc: |
  Converts GTEx rMATS splice event output to mean junction count matrices by GTEx subgroup

requirements:
  DockerRequirement:
    dockerPull: "pgc-images.sbgenomics.com/childrens-bti/impact_trial_cwl:latest"
  InitialWorkDirRequirement:
    listing:
      - $(inputs.gtex_se_rmats)
      - $(inputs.gtex_ri_rmats)
      - $(inputs.gtex_a3ss_rmats)
      - $(inputs.gtex_a5ss_rmats)
      - $(inputs.histologies_file)
      - $(inputs.gtex_read_counts)

baseCommand: [Rscript, --vanilla, /opt/impact_trial/scripts/create-gtex-matrices/create-gtex-jc-count-matrices.R]

inputs:
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
  
  gtex_read_counts:
    type: File
    doc: "GTEx input read counts file"

outputs:
  gtex_se_junction_matrix:
    type: File
    outputBinding:
      glob: "results/gtex-se-junction-mat-raw.qs2"
  
  gtex_se_target_matrix:
    type: File
    outputBinding:
      glob: "results/gtex-se-target-mat-raw.qs2"
  
  gtex_ri_junction_matrix:
    type: File
    outputBinding:
      glob: "results/gtex-ri-junction-mat-raw.qs2"
  
  gtex_ri_target_matrix:
    type: File
    outputBinding:
      glob: "results/gtex-ri-target-mat-raw.qs2"
  
  gtex_a3ss_junction_matrix:
    type: File
    outputBinding:
      glob: "results/gtex-a3ss-junction-mat-raw.qs2"
  
  gtex_a5ss_junction_matrix:
    type: File
    outputBinding:
      glob: "results/gtex-a5ss-junction-mat-raw.qs2"
