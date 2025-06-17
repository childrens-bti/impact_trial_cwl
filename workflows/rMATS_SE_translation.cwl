class: "CommandLineTool"
cwlVersion: "v1.0"
baseCommand: ["Rscript"]
inputs:
  gtf_file:
    type: "File"
    inputBinding:
      position: 1
      prefix: "--gtf"
  genome_fa:
    type: "File"
    inputBinding:
      position: 2
      prefix: "--genome"
  rMATS_file:
    type: "File"
    inputBinding:
      position: 3
      prefix: "--rmats"
  output_fasta:
    type: "File"
    outputBinding:
      glob: "rMATS_SE_peptides_test.fasta"
outputs:
  peptides:
    type: "File"
    outputBinding:
      glob: "rMATS_SE_peptides_test.fasta"
requirements:
  - class: "InlineJavascriptRequirement"
  - class: "DockerRequirement"
    dockerPull: "rocker/r-ver:4.1.0"
    env:
      R_LIBS_USER: "/usr/local/lib/R/site-library"
      R_ENVIRON: "/usr/local/lib/R/etc/Renviron"