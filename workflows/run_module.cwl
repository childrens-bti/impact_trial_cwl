class: CommandLineTool
cwlVersion: v1.0
baseCommand: ["bash", "scripts/translation/run_module.sh"]
inputs:
  snv_input:
    type: File
    inputBinding:
      position: 1
      prefix: "-i"
  snv_output:
    type: File
    outputBinding:
      glob: "analyses/translation/results/snv-consensus-plus-hotspots.maf.fasta"
  star_fusion_input:
    type: File
    inputBinding:
      position: 2
      prefix: "-i"
  star_fusion_output:
    type: File
    outputBinding:
      glob: "analyses/translation/results/fusion-starfusion-peptides.fasta"
  arriba_fusion_input:
    type: File
    inputBinding:
      position: 3
      prefix: "-i"
  arriba_fusion_output:
    type: File
    outputBinding:
      glob: "analyses/translation/results/fusion-arriba-peptides.fasta"
outputs:
  snv_output:
    type: File
    outputBinding:
      glob: "analyses/translation/results/snv-consensus-plus-hotspots.maf.fasta"
  star_fusion_output:
    type: File
    outputBinding:
      glob: "analyses/translation/results/fusion-starfusion-peptides.fasta"
  arriba_fusion_output:
    type: File
    outputBinding:
      glob: "analyses/translation/results/fusion-arriba-peptides.fasta"