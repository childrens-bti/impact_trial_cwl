id: translation
label: Translation Analysis
class: Workflow
inputs:
  snv_input:
    type: File
    inputBinding:
      position: 1
  star_fusion_input:
    type: File
    inputBinding:
      position: 2
  arriba_fusion_input:
    type: File
    inputBinding:
      position: 3
  digestion:
    type: string
    default: "lysC"
    inputBinding:
      position: 4
outputs:
  snv_output:
    type: File
    outputBinding:
      glob: "$(inputs.snv_input.basename).fasta"
  star_fusion_output:
    type: File
    outputBinding:
      glob: "$(inputs.star_fusion_input.basename).fasta"
  arriba_fusion_output:
    type: File
    outputBinding:
      glob: "$(inputs.arriba_fusion_input.basename).fasta"

steps:
  snv_translation:
    run: 01-snv-translate.cwl
    in:
      input: snv_input
      digestion: digestion
    out: [snv_output]

  star_fusion_translation:
    run: 02-star-fusion-translate.cwl
    in:
      input: star_fusion_input
      digestion: digestion
    out: [star_fusion_output]

  arriba_fusion_translation:
    run: 03-arriba-fusion-translate.cwl
    in:
      input: arriba_fusion_input
      digestion: digestion
    out: [arriba_fusion_output]