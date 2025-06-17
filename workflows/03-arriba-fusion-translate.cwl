id: 03-arriba-fusion-translate
label: Arriba Fusion Translate
class: CommandLineTool
baseCommand: python
inputs:
  - id: input_file
    type: File
    inputBinding:
      position: 1
      prefix: -i
  - id: output_file
    type: File
    inputBinding:
      position: 2
      prefix: -o
  - id: digestion_enzyme
    type: string
    default: lysC
    inputBinding:
      prefix: -d
outputs:
  - id: output_peptides
    type: File
    outputBinding:
      glob: $(inputs.output_file)
requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: pgc-images.sbgenomics.com/rokita-lab/impact-trial:latest
    dockerImageId: pgc-images.sbgenomics.com/rokita-lab/impact-trial:latest
    dockerOutput: true
    dockerCommand: python
    dockerArgs:
      - -i
      - $(inputs.input_file)
      - -o
      - $(inputs.output_file)
      - -d
      - $(inputs.digestion_enzyme)