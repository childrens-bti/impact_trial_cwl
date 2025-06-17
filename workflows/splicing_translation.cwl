id: splicing_translation
label: Splicing Translation
class: Workflow
inputs:
  gtf_file:
    type: File
    inputBinding:
      position: 1
  genome_fa:
    type: File
    inputBinding:
      position: 2
  rMATS_file:
    type: File
    inputBinding:
      position: 3
  output_fasta:
    type: File
    inputBinding:
      position: 4
steps:
  rMATS_SE_translation:
    run: rMATS_SE_translation.cwl
    in:
      gtf_file: gtf_file
      genome_fa: genome_fa
      rMATS_file: rMATS_file
      output_fasta: output_fasta
    out: [output_fasta]