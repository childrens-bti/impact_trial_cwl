# Tumor-specific splice product translation
#
# Ryan Corbett / Bicna Song
#
# July 2025

# This script performs in silico translation of tumor-specific splice products in an ORF-unaware manner to retain the maximum number of potentially translated peptides that meet a specified AA length threshold. 

### Load packages
library(rprojroot)
library(tidyverse)
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)
library(glue)
library(optparse)

### Parse command-line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "tumor-specific splice event file", metavar = "input"),
  make_option(c("--min_peptide_length"), type = "integer", default = 8, help = "Minimum peptide length [default %default]", metavar = "length"),
  make_option(c("--n_flank_bases"), type = "integer", default = 20, help = "Number of junction-flanking nucleotides to retain for translation [default %default]", metavar = "n_flank_bases"),
  make_option(c("--sample_id"), type = "character", default = NULL, help = "Sample ID to append to output file", metavar = "sample_id")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser, args = commandArgs(trailingOnly = TRUE))
event_file <- opt$input
min_peptide_length <- opt$min_peptide_length
n_flank_bases <- opt$n_flank_bases
sample <- opt$sample_id

### Set directory paths
root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "translation")
results_dir <- file.path(analysis_dir, "results")

source(file.path(analysis_dir, "util", "get_seq.R"))
source(file.path(analysis_dir, "util", "translate_all_phases.R"))

### Input genome and output
genome_fa <- file.path(data_dir, "GRCh38.primary_assembly.genome.fa.gz")
output_fasta <- file.path(results_dir, glue::glue("{sample}_tumor_specific_splice_event_peptides.fasta"))

genome <- readDNAStringSet(genome_fa)
names(genome) <- sub(" .*$", "", names(genome))

### Load files
events <- read_tsv(event_file, col_types = cols()) 

# filter for splice event junctions
junctions <- events %>%
  dplyr::filter(type == "junction") %>%
  dplyr::mutate(
    # split junction coordinates
    chr = unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[1])),
    up_exon_start = as.numeric(unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[2]))),
    jc_up_coord = as.numeric(unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[3]))),
    jc_down_coord = as.numeric(unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[4]))),
    down_exon_end = as.numeric(unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[5])))
  ) %>%
  dplyr::mutate(up_exon_start = up_exon_start + 1,
                jc_down_coord = jc_down_coord + 1)

set.seed(123)

### Translate junctions to peptides 

peptide_list <- list()

junction_report_log <- tibble(event = junctions$coordinates,
                     type = junctions$type,
                     strand = junctions$strand,
                     n_jc_peptides = rep(0, nrow(junctions)))

for (i in 1:nrow(junctions)) {
  event <- junctions[i, ]
  chr <- as.character(event$chr)
  strand <- as.character(event$strand)
  coord_str <- event$coordinates
  gene <- as.character(event$geneSymbol)
  
  # obtain junction intervals
  
  intervals <- list(
    up_jc = c("start" = ifelse(
      n_flank_bases > event$jc_up_coord - event$up_exon_start,
      event$up_exon_start,
      event$jc_up_coord - n_flank_bases),
      "end" = as.numeric(event$jc_up_coord)),
    down_jc = c("start" = as.numeric(event$jc_down_coord),
                "end" = ifelse(
                  n_flank_bases > event$down_exon_end - event$jc_down_coord,
                  event$down_exon_end,
                  event$jc_down_coord + n_flank_bases))
  )
  
  # obtain nucleotide sequences 
  seqs <- map(intervals, ~ get_seq(event$chr, .x["start"], .x["end"], event$strand))
  
  # get junction sequence
    if (strand == "+"){
      
      junction_seq <- xscat(seqs$up_jc,
                             seqs$down_jc)
      
    } else if (strand == "-"){
      
      junction_seq <- xscat(seqs$down_jc,
                               seqs$up_jc)
      
    }
  
  # define variable to tally number of incl and skip peptides translated and passing filter
  peptide_ct = 0
  
  # translate in three phases
  peptides <- translate_all_phases(junction_seq,
                                   coord_str = coord_str,
                                   gene = gene)
  
  # filter for min peptide length
  filtered_peptides <- peptides %>%
    map(~ .x[width(.x) >= min_peptide_length]) %>%
    keep(~ length(.x) > 0)
  
  # add peptides to list
  peptide_list <- c(peptide_list, filtered_peptides)
  
  # get number of peptides for jc 
  peptide_ct = length(filtered_peptides)

  # add peptide n to report log
  junction_report_log$n_jc_peptides[i] <- peptide_ct
  
} 


# filter for splice event targets
targets <- events %>%
  dplyr::filter(type != "junction") %>%
  dplyr::mutate(
    # split junction coordinates
    chr = unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[1])),
    target_start = as.numeric(unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[2]))),
    target_end = as.numeric(unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[3])))
  ) %>%
  dplyr::mutate(target_start = target_start + 1,
                target_end = target_end + 1)


### Translate targets to peptides 
target_report_log <- tibble(event = targets$coordinates,
                              type = targets$type,
                              strand = targets$strand,
                              n_tgt_peptides = rep(0, nrow(targets)))

for (i in 1:nrow(targets)) {
  event <- targets[i, ]
  chr <- as.character(event$chr)
  strand <- as.character(event$strand)
  coord_str <- event$coordinates
  gene <- as.character(event$geneSymbol)
  
  # obtain junction intervals
  
  intervals <- list(
    target = c("start" = event$target_start,
               "end" = event$target_end)
  )
  
  # obtain nucleotide sequences 
  seq <- map(intervals, ~ get_seq(event$chr, .x["start"], .x["end"], event$strand))
  
  # get junction sequence
  if (strand == "+"){

    target_seq <- xscat(seq$target)

  } else if (strand == "-"){

    target_seq <- xscat(seq$target)

  }
  
  # define variable to tally number of incl and skip peptides translated and passing filter
  peptide_ct = 0
  
  # translate in three phases
  peptides <- translate_all_phases(target_seq,
                                   coord_str = coord_str,
                                   gene = gene)
  
  # filter for min peptide length
  filtered_peptides <- peptides %>%
    map(~ .x[width(.x) >= min_peptide_length]) %>%
    keep(~ length(.x) > 0)
  
  # add peptides to list
  peptide_list <- c(peptide_list, filtered_peptides)
  
  # get number of peptides for jc 
  peptide_ct = length(filtered_peptides)
  
  # add peptide n to report log
  target_report_log$n_tgt_peptides[i] <- peptide_ct
  
} 

# Merge junction- and target-derived peptides, and write to output
peptides <- do.call(c, peptide_list)
writeXStringSet(peptides, filepath = output_fasta, format = "fasta")

write_tsv(junction_report_log,
          file.path(results_dir, glue::glue("{sample}-tumor-specific-splice-junction-translation-report.tsv")))
write_tsv(junction_report_log,
          file.path(results_dir, glue::glue("{sample}-tumor-specific-splice-target-translation-report.tsv")))

# Session info

sessionInfo()

