
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

baseCommand: [bash, -c]

arguments:
  - valueFrom: |
      set -e
      # Build BLAST database if needed
      if [ ! -f "uniprot_db.pin" ]; then
          makeblastdb -in $(inputs.uniprot_fasta.path) -dbtype prot -out uniprot_db
      fi
      
      # Run BLAST
      blastp -query $(inputs.custom_fasta.path) \
             -db uniprot_db \
             -out blast_results.tsv \
             -outfmt "6 qseqid sseqid pident length qlen slen"
      
      # Extract 100% matches
      awk '$3 == 100 && $4 == $5 {print $1}' blast_results.tsv | sort -u > full_matches.txt
      
      # Count statistics
      total=$(grep -c "^>" $(inputs.custom_fasta.path))
      matches=$(wc -l < full_matches.txt)
      remaining=$(( total - matches ))
      
      echo "Total peptides: $total"
      echo "Filtered out (100% UniProt matches): $matches"
      echo "Remaining peptides: $remaining"
      
      # Filter peptides using seqkit
      seqkit grep -v -f full_matches.txt $(inputs.custom_fasta.path) > custom_filtered.fasta

inputs:
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
      glob: "custom_filtered.fasta"
  
  blast_results:
    type: File
    outputBinding:
      glob: "blast_results.tsv"
