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
# Use current working directory (CWL stages files here)
results_dir <- file.path("results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Source utility functions from Docker image
source("/opt/impact_trial/scripts/translation/util/get_seq.R")
source("/opt/impact_trial/scripts/translation/util/translate_all_phases.R")

### Input genome and output (file staged in current directory by CWL)
genome_fa <- "GRCh38.primary_assembly.genome.fa.gz"
output_fasta <- file.path(results_dir, glue::glue("{sample}_tumor_specific_splice_event_peptides.fasta"))

genome <- readDNAStringSet(genome_fa)
names(genome) <- sub(" .*$", "", names(genome))

### Load files
events <- read_tsv(event_file, col_types = cols()) %>%
  dplyr::mutate(
    chr = unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[1])),
    jc_up_coord = as.numeric(unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[3]))),
    jc_down_coord = as.numeric(unlist(lapply(str_split(coordinates, ":|-|_"), function(x) x[4])))
  ) %>%
  dplyr::mutate(jc_down_coord = jc_down_coord + 1)

set.seed(123)

### Translate peptides (core loop)

peptide_list <- list()

report_log <- tibble(event = events$coordinates,
                     type = events$type,
                     strand = events$strand,
                     n_jc_peptides = rep(0, nrow(events)))

for (i in 1:nrow(events)) {
  event <- events[i, ]
  chr <- as.character(event$chr)
  strand <- as.character(event$strand)
  coord_str <- event$coordinates
  gene <- as.character(event$geneSymbol)
  
  # obtain junction intervals
  
  intervals <- list(
    up_jc = c("start" = as.numeric(event$jc_up_coord - n_flank_bases),
                "end" = as.numeric(event$jc_up_coord)),
    down_jc = c("start" = as.numeric(event$jc_down_coord),
                "end" = as.numeric(event$jc_down_coord + n_flank_bases))
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
  report_log$n_jc_peptides[i] <- peptide_ct
  
} 

# Write peptides to output
peptides <- do.call(c, peptide_list)
writeXStringSet(peptides, filepath = output_fasta, format = "fasta")

write_tsv(report_log,
          file.path(results_dir, glue::glue("{sample}-tumor-specific-splice-product-translation-report.tsv")))

# Session info

sessionInfo()

