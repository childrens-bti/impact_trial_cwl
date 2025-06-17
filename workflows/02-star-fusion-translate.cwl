id: star_fusion_translate
label: STAR Fusion Translate
class: CommandLineTool
baseCommand: python
inputs:
  - id: input_file
    type: File
    inputBinding:
      position: 1
  - id: output_file
    type: File
    inputBinding:
      position: 2
  - id: digestion_enzyme
    type: string
    default: lysC
    inputBinding:
      position: 3
outputs:
  - id: output_peptides
    type: File
    outputBinding:
      glob: $(inputs.output_file)