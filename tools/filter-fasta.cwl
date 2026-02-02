cwlVersion: v1.2
class: CommandLineTool

label: Filter FASTA against UniProt
doc: |
  BLASTs custom peptides against UniProt canonical database and filters out 100% matches

requirements:
  DockerRequirement:
    dockerPull: "pgc-images.sbgenomics.com/childrens-bti/impact_trial_cwl:latest"
  InitialWorkDirRequirement:
    listing:
      - $(inputs.uniprot_fasta)
      - $(inputs.custom_fasta)

baseCommand: [bash, /opt/impact_trial/scripts/translation/06-filter-fasta.sh]

arguments:
  - prefix: --custom_fasta
    valueFrom: $(inputs.custom_fasta.path)
  - prefix: --uniprot_fasta
    valueFrom: $(inputs.uniprot_fasta.path)
  - prefix: --dna_sample_id
    valueFrom: $(inputs.dna_sample_id)
  - prefix: --rna_sample_id
    valueFrom: $(inputs.rna_sample_id)

inputs:
  dna_sample_id:
    type: string
    doc: "DNA sample identifier"
  
  rna_sample_id:
    type: string
    doc: "RNA sample identifier"
  
  custom_fasta:
    type: File
    doc: "Merged custom peptide FASTA"
  
  uniprot_fasta:
    type: File
    doc: "UniProt canonical proteome FASTA"

outputs:
  filtered_fasta:
    type: File
    outputBinding:
      glob: "$(inputs.dna_sample_id)_$(inputs.rna_sample_id)_tumor_specific_filtered.fasta"
  
  blast_results:
    type: File
    outputBinding:
      glob: "$(inputs.dna_sample_id)_$(inputs.rna_sample_id)_uniprot_blast_results.tsv"
  
  full_matches:
    type: File
    outputBinding:
      glob: "$(inputs.dna_sample_id)_$(inputs.rna_sample_id)_uniprot_full_matches.txt"
