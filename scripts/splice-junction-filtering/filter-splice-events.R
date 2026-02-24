## Filter splice events for specificity relative to GTEx
##
## Ryan Corbett
##
## Nov 2025

# This script loads tumors rMATS files, pull splice event junction and target coordinates, identifies those that are tumor-specific by comparing to GTEx splice events junctions/targets

library(tidyverse)
library(qs2)
library(optparse)

### Parse command-line options
option_list <- list(
  make_option(c("--tumor_se_rmats"), type = "character", default = NULL, help = "Tumor rMATS SE JC file", metavar = "tumor_se_rmats"),
  make_option(c("--tumor_ri_rmats"), type = "character", default = NULL, help = "Tumor rMATS RI JC file", metavar = "tumor_ri_rmats"),
  make_option(c("--tumor_a3ss_rmats"), type = "character", default = NULL, help = "Tumor rMATS A3SS JC file", metavar = "tumor_a3ss_rmats"),
  make_option(c("--tumor_a5ss_rmats"), type = "character", default = NULL, help = "Tumor rMATS A5SS JC file", metavar = "tumor_a5ss_rmats"),
  make_option(c("--sample_id"), type = "character", default = NULL, help = "Sample ID", metavar = "sample_id"),
  make_option(c("--se_gtex_junction_mat_file"), type = "character", default = NULL, help = "GTEx SE junction matrix", metavar = "se_gtex_junction_mat_file"),
  make_option(c("--se_gtex_target_mat_file"), type = "character", default = NULL, help = "GTEx SE target matrix", metavar = "se_gtex_target_mat_file"),
  make_option(c("--ri_gtex_junction_mat_file"), type = "character", default = NULL, help = "GTEx RI junction matrix", metavar = "ri_gtex_junction_mat_file"),
  make_option(c("--ri_gtex_intron_mat_file"), type = "character", default = NULL, help = "GTEx RI target matrix", metavar = "ri_gtex_intron_mat_file"),
  make_option(c("--a3ss_gtex_junction_mat_file"), type = "character", default = NULL, help = "GTEx A3SS junction matrix", metavar = "a3ss_gtex_junction_mat_file"),
  make_option(c("--a5ss_gtex_junction_mat_file"), type = "character", default = NULL, help = "GTEx A5SS junction matrix", metavar = "a5ss_gtex_junction_mat_file"),
  make_option(c("--gencode_gtf"), type = "character", default = NULL, help = "GENCODE GTF", metavar = "gencode_gtf")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser, args = commandArgs(trailingOnly = TRUE))
tumor_se_jc_file <- opt$tumor_se_rmats
tumor_ri_jc_file <- opt$tumor_ri_rmats
tumor_a3ss_jc_file <- opt$tumor_a3ss_rmats
tumor_a5ss_jc_file <- opt$tumor_a5ss_rmats
sample_id <- opt$sample_id
se_gtex_junction_mat_file <- opt$se_gtex_junction_mat_file
se_gtex_target_mat_file <- opt$se_gtex_target_mat_file
ri_gtex_junction_mat_file <- opt$ri_gtex_junction_mat_file
ri_gtex_intron_mat_file <- opt$ri_gtex_intron_mat_file
a3ss_gtex_junction_mat_file <- opt$a3ss_gtex_junction_mat_file
a5ss_gtex_junction_mat_file <- opt$a5ss_gtex_junction_mat_file
gencode_gtf_file <- opt$gencode_gtf

required_args <- c(
  "tumor_se_rmats" = tumor_se_jc_file,
  "tumor_ri_rmats" = tumor_ri_jc_file,
  "tumor_a3ss_rmats" = tumor_a3ss_jc_file,
  "tumor_a5ss_rmats" = tumor_a5ss_jc_file,
  "sample_id" = sample_id,
  "se_gtex_junction_mat_file" = se_gtex_junction_mat_file,
  "se_gtex_target_mat_file" = se_gtex_target_mat_file,
  "ri_gtex_junction_mat_file" = ri_gtex_junction_mat_file,
  "ri_gtex_intron_mat_file" = ri_gtex_intron_mat_file,
  "a3ss_gtex_junction_mat_file" = a3ss_gtex_junction_mat_file,
  "a5ss_gtex_junction_mat_file" = a5ss_gtex_junction_mat_file,
  "gencode_gtf" = gencode_gtf_file
)

missing_args <- names(required_args)[vapply(required_args, function(x) is.null(x) || x == "", logical(1))]
if (length(missing_args) > 0) {
  stop(glue::glue("Missing required arguments: {paste(missing_args, collapse = ', ')}"))
}

### Set directory paths
# Input/output paths are provided by CWL arguments and resolved by the runtime.

## functions for determining splice event specificity

# tumor-specific junctions
get_tumor_specific_junctions <- function(df,
                                         gtex_mat){

  spec_jc_df <- df %>%
    dplyr::mutate(boundary = unlist(lapply(str_split(junction, "-"), function(x) x[2]))) %>%
    # join gtex junction counts
    left_join(gtex_mat) %>%
    dplyr::mutate(
      # N gtex subgroups w/ mean count < 1 (exclude testis)
      n_low_ct_gtex = rowSums(across(contains(c("Brain", "Nerve", "Skin", "Testis")), ~ . < 1, .names = NULL), na.rm = TRUE),
      # N subgroups NA
      n_na_gtex = rowSums(across(contains(c("Brain", "Nerve", "Skin", "Testis")), ~ is.na(.), .names = NULL), na.rm = TRUE)) %>%
    
    # determine if junction is tumor-specific
    dplyr::mutate(junction_preference = case_when(
      # Low ct in all gtex groups with non-NA values
      # Tumor junction ct >= 10
      n_low_ct_gtex == n_gtex_groups - n_na_gtex & junction_count >= 10 ~ "Tumor-specific junction",
      TRUE ~ "No specificity"
    ))

  # pull all boundaries called non-specific 
  nonspecific_boundaries <- spec_jc_df %>%
    dplyr::filter(junction_preference == "No specificity") %>%
    pull(boundary)
  
  # retain only tumor-specific junctions
  spec_jc_df <- spec_jc_df %>%
    dplyr::filter(junction_preference == "Tumor-specific junction",
                  # remove cases where boundary called nonspecific as part of another junction
                  !boundary %in% nonspecific_boundaries)
  
  return(spec_jc_df)
  
}

# tumor-specific targets 
get_tumor_specific_targets <- function(df,
                                       gtex_mat){

  spec_target_df <- df %>%
    # join gtex junction counts
    left_join(gtex_mat) %>%
    dplyr::mutate(
      # N gtex subgroups w/ mean count < 1
      n_low_ct_gtex = rowSums(across(contains(c("Brain", "Nerve", "Skin", "Testis")), ~ . < 1, .names = NULL), na.rm = TRUE),
      # N subgroups NA
      n_na_gtex = rowSums(across(contains(c("Brain", "Nerve", "Skin", "Testis")), ~ is.na(.), .names = NULL), na.rm = TRUE)) %>%
    
    # determine if target is tumor-specific
    dplyr::mutate(target_preference = case_when(
      # Low ct in all gtex groups with non-NA values
      # Tumor junction ct >= 10
      n_low_ct_gtex == n_gtex_groups - n_na_gtex & target_count >= 10 ~ "Tumor-specific target",
      TRUE ~ "No specificity"
    )) %>%
    # retain only tumor-specific targets
    dplyr::filter(target_preference == "Tumor-specific target")
  
  return(spec_target_df)
  
}

# Load SE rMATS file and define splice, target, and junction IDs
tumor_se_rmats <- read_tsv(tumor_se_jc_file) %>%
  dplyr::mutate(up_incl_jc = glue::glue("{chr}:{upstreamES}-{upstreamEE}_{exonStart_0base}-{exonEnd}"),
                down_incl_jc = glue::glue("{chr}:{exonStart_0base}-{exonEnd}_{downstreamES}-{downstreamEE}"),
                target = glue::glue("{chr}:{exonStart_0base}_{exonEnd}"),
                skip_jc = glue::glue("{chr}:{upstreamES}-{upstreamEE}_{downstreamES}-{downstreamEE}"),
                sample_id = sample_id) %>%
  # select relevant columns
  dplyr::select(sample_id, geneSymbol, strand,
                up_incl_jc,
                down_incl_jc, target, skip_jc,
                upstream_to_target_count,
                target_to_downstream_count,
                target_count,
                upstream_to_downstream_count) %>%
  dplyr::filter(upstream_to_target_count != 0 | target_to_downstream_count != 0 | upstream_to_downstream_count != 0 | target != 0)

# Load SE GTEx junction & target matrices 
gtex_junction_df <- qs2::qs_read(se_gtex_junction_mat_file)

gtex_target_df <- qs2::qs_read(se_gtex_target_mat_file)

# define number of GTEx groups (exclude Testis)
n_gtex_groups <- sum(grepl("Brain|Nerve|Skin|Testis", names(gtex_junction_df)))

# create data frame of all junctions
se_junction_mat <- tumor_se_rmats %>%
  # filter for up_incl junctions
  dplyr::select(sample_id, geneSymbol, strand,
                up_incl_jc, upstream_to_target_count) %>%
  dplyr::rename(junction = up_incl_jc,
                junction_ct = upstream_to_target_count) %>%
  # add down_incl junctions
  bind_rows(tumor_se_rmats %>%
              dplyr::select(sample_id, geneSymbol, strand,
                            down_incl_jc, target_to_downstream_count) %>%
              dplyr::rename(junction = down_incl_jc,
                            junction_ct = target_to_downstream_count)) %>%
  # add skip junctions
  bind_rows(tumor_se_rmats %>%
              dplyr::select(sample_id, geneSymbol, strand,
                            skip_jc, upstream_to_downstream_count) %>%
              dplyr::rename(junction = skip_jc,
                            junction_ct = upstream_to_downstream_count)) %>%
  distinct(junction, .keep_all = TRUE) %>%
  # for rare cases of duplicated junctions, calculate median count
  group_by(sample_id, geneSymbol, strand, junction) %>%
  summarise(junction_count = median(junction_ct)) %>%
  ungroup() %>%
  # only retain junctions with non-zero counts
  dplyr::filter(junction_count != 0)

# calculate junction specificity
se_junction_df <- get_tumor_specific_junctions(se_junction_mat,
                                               gtex_junction_df)

# create data frame of all targets (here, exons)
se_target_mat <- tumor_se_rmats %>%
  # filter for up_incl junctions
  dplyr::select(sample_id, geneSymbol, strand,
                target, target_count) %>%
  # for rare cases of duplicated junctions, calculate median count
  group_by(sample_id, geneSymbol, strand, target) %>%
  summarise(target_count = median(target_count)) %>%
  ungroup()

# calculate junction specificity
se_target_df <- get_tumor_specific_targets(se_target_mat,
                                           gtex_target_df)

# Load RI rMATS file and define splice, target, and junction IDs
tumor_ri_rmats <- read_tsv(tumor_ri_jc_file) %>%
  dplyr::mutate(up_incl_jc = glue::glue("{chr}:{upstreamES}-{upstreamEE}_{upstreamEE}-{downstreamES}"),
                down_incl_jc = glue::glue("{chr}:{upstreamEE}-{downstreamES}_{downstreamES}-{downstreamEE}"),
                target = glue::glue("{chr}:{upstreamEE}_{downstreamES}"),
                skip_jc = glue::glue("{chr}:{upstreamES}-{upstreamEE}_{downstreamES}-{downstreamEE}"),
                sample_id = sample_id) %>%
  dplyr::rename(target_count = intron_count) %>%
  # select relevant columns
  dplyr::select(sample_id, geneSymbol, strand,
                up_incl_jc,
                down_incl_jc, target, skip_jc,
                upstream_to_intron_count,
                intron_to_downstream_count,
                target_count,
                upstream_to_downstream_count) %>%
  dplyr::filter(upstream_to_intron_count != 0 | intron_to_downstream_count != 0 | upstream_to_downstream_count != 0 | target != 0)

# Load GTEx data 
ri_gtex_junction_df <- qs2::qs_read(ri_gtex_junction_mat_file)

gtex_intron_df <- qs2::qs_read(ri_gtex_intron_mat_file)

# create data frame of all junctions in rMATS
ri_junction_mat <- tumor_ri_rmats %>%
  # filter for up_incl junctions
  dplyr::select(sample_id, geneSymbol, strand,
                up_incl_jc, upstream_to_intron_count) %>%
  dplyr::rename(junction = up_incl_jc,
                junction_ct = upstream_to_intron_count) %>%
  # add down_incl junctions
  bind_rows(tumor_ri_rmats %>%
              dplyr::select(sample_id, geneSymbol, strand,
                            down_incl_jc, intron_to_downstream_count) %>%
              dplyr::rename(junction = down_incl_jc,
                            junction_ct = intron_to_downstream_count)) %>%
  # add skip junctions
  bind_rows(tumor_ri_rmats %>%
              dplyr::select(sample_id, geneSymbol, strand,
                            skip_jc, upstream_to_downstream_count) %>%
              dplyr::rename(junction = skip_jc,
                            junction_ct = upstream_to_downstream_count)) %>%
  distinct(junction, .keep_all = TRUE) %>%
  # for rare cases of duplicated junctions, calculate median count
  group_by(sample_id, geneSymbol, strand, junction) %>%
  summarise(junction_count = median(junction_ct)) %>%
  ungroup() %>%
  # remove junctions with zero counts
  dplyr::filter(!junction %in% se_junction_mat$junction,
                junction_count > 0)

# calculate junction specificity
ri_junction_df <- get_tumor_specific_junctions(ri_junction_mat,
                                               ri_gtex_junction_df)

# create data frame of all targets (here, introns) in rMATS
ri_target_mat <- tumor_ri_rmats %>%
  dplyr::select(sample_id, geneSymbol, strand,
                target, target_count) %>%
  distinct(target, .keep_all = TRUE) %>%
  # for rare cases of duplicated targets, calculate median count
  group_by(sample_id, geneSymbol, strand,
           target) %>%
  summarise(target_count = median(target_count)) %>%
  ungroup() %>% 
  # retain introns with non-zero counts
  dplyr::filter(target_count > 0)

# calculate target specificity
ri_target_df <- get_tumor_specific_targets(ri_target_mat,
                                           gtex_intron_df)

# A3SS events
tumor_a3ss_rmats <- read_tsv(tumor_a3ss_jc_file) %>%
  dplyr::filter(long_to_flanking_count != 0 | short_to_flanking_count != 0) %>%
  dplyr::mutate(long_incl_jc = case_when(
    strand == "+" ~ glue::glue("{chr}:{flankingES}-{flankingEE}_{longExonStart_0base}-{longExonEnd}"),
    strand == "-" ~ glue::glue("{chr}:{longExonStart_0base}-{longExonEnd}_{flankingES}-{flankingEE}"))) %>%
  dplyr::mutate(short_incl_jc = case_when(
    strand == "+" ~ glue::glue("{chr}:{flankingES}-{flankingEE}_{shortES}-{shortEE}"),
    strand == "-" ~ glue::glue("{chr}:{shortES}-{shortEE}_{flankingES}-{flankingEE}")
  )) %>%
  dplyr::mutate(sample_id = sample_id) %>%
  # select relevant columns
  dplyr::select(sample_id, geneSymbol, strand,
                long_incl_jc,
                short_incl_jc,
                long_to_flanking_count,
                short_to_flanking_count)

# Load GTEx data 
a3ss_gtex_junction_df <- qs2::qs_read(a3ss_gtex_junction_mat_file)

# create df of merged junctions 
a3ss_junction_mat <- tumor_a3ss_rmats %>%
  # filter for long incl junctions
  dplyr::select(sample_id, geneSymbol, strand,
                long_incl_jc, long_to_flanking_count) %>%
  dplyr::rename(junction = long_incl_jc,
                junction_ct = long_to_flanking_count) %>%
  # add short incl junctions
  bind_rows(tumor_a3ss_rmats %>%
              dplyr::select(sample_id, geneSymbol, strand,
                            short_incl_jc, short_to_flanking_count) %>%
              dplyr::rename(junction = short_incl_jc,
                            junction_ct = short_to_flanking_count)) %>%
  distinct(junction, .keep_all = TRUE) %>%
  # for rare cases of duplicated junctions, calculate median count
  group_by(sample_id, geneSymbol, strand, junction) %>%
  summarise(junction_count = median(junction_ct)) %>%
  ungroup() %>%
  dplyr::filter(!junction %in% c(se_junction_mat$junction,
                                 ri_junction_mat$junction),
                junction_count > 0)

# calculate junction specificity
a3ss_junction_df <- get_tumor_specific_junctions(a3ss_junction_mat,
                                                 a3ss_gtex_junction_df)

# A5SS
tumor_a5ss_rmats <- read_tsv(tumor_a5ss_jc_file) %>%
  dplyr::filter(long_to_flanking_count != 0 | short_to_flanking_count != 0) %>%
  dplyr::mutate(long_incl_jc = case_when(
    strand == "+" ~ glue::glue("{chr}:{longExonStart_0base}-{longExonEnd}_{flankingES}-{flankingEE}"),
    strand == "-" ~ glue::glue("{chr}:{flankingES}-{flankingEE}_{longExonStart_0base}-{longExonEnd}"))) %>%
  dplyr::mutate(short_incl_jc = case_when(
    strand == "+" ~ glue::glue("{chr}:{shortES}-{shortEE}_{flankingES}-{flankingEE}"),
    strand == "-" ~ glue::glue("{chr}:{flankingES}-{flankingEE}_{shortES}-{shortEE}")
  )) %>%
  dplyr::mutate(sample_id = sample_id) %>%
  # select relevant columns
  dplyr::select(sample_id, geneSymbol, strand,
                long_incl_jc,
                short_incl_jc,
                long_to_flanking_count,
                short_to_flanking_count)

# Load GTEx data 
a5ss_gtex_junction_df <- qs2::qs_read(a5ss_gtex_junction_mat_file)

# create df of merged junctions
a5ss_junction_mat <- tumor_a5ss_rmats %>%
  # filter for long incl junctions
  dplyr::select(sample_id, geneSymbol, strand,
                long_incl_jc, long_to_flanking_count) %>%
  dplyr::rename(junction = long_incl_jc,
                junction_ct = long_to_flanking_count) %>%
  # add short incl junctions
  bind_rows(tumor_a5ss_rmats %>%
              dplyr::select(sample_id, geneSymbol, strand,
                            short_incl_jc, short_to_flanking_count) %>%
              dplyr::rename(junction = short_incl_jc,
                            junction_ct = short_to_flanking_count)) %>%
  distinct(junction, .keep_all = TRUE) %>%
  # for rare cases of duplicated junctions, calculate median count
  group_by(sample_id, geneSymbol, strand, junction) %>%
  summarise(junction_count = median(junction_ct)) %>%
  ungroup() %>%
  dplyr::filter(!junction %in% c(se_junction_mat$junction,
                                 ri_junction_mat$junction,
                                 a3ss_junction_mat$junction),
                junction_count > 0)

# calculate junction specificity
a5ss_junction_df <- get_tumor_specific_junctions(a5ss_junction_mat,
                                                 a5ss_gtex_junction_df)

# Define splice event preferences

# A3 loss events defined by a3 short inclusion junction 
a3_loss_events <- tumor_a3ss_rmats %>%
  dplyr::filter(short_to_flanking_count >= 5) %>%
  dplyr::select(short_incl_jc) %>%
  dplyr::rename(junction = short_incl_jc) %>%
  distinct() %>%
  dplyr::mutate(type = "A3-") 

# A3 gain events defined by a3 long inclusion junction
a3_gain_events <- tumor_a3ss_rmats %>%
  dplyr::filter(long_to_flanking_count >= 5) %>%
  dplyr::select(long_incl_jc) %>%
  dplyr::rename(junction = long_incl_jc) %>%
  distinct() %>%
  dplyr::mutate(type = "A3+") 

# A5 loss events defined by a5 short inclusion junction
a5_loss_events <- tumor_a5ss_rmats %>%
  dplyr::filter(short_to_flanking_count >= 5) %>%
  dplyr::select(short_incl_jc) %>%
  dplyr::rename(junction = short_incl_jc) %>%
  distinct() %>%
  dplyr::mutate(type = "A5-") 

# A5 gain events defined by a5 long inclusion junction
a5_gain_events <- tumor_a5ss_rmats %>%
  dplyr::filter(long_to_flanking_count >= 5) %>%
  dplyr::select(long_incl_jc) %>%
  dplyr::rename(junction = long_incl_jc) %>%
  distinct() %>%
  dplyr::mutate(type = "A5+") 

# define exon inclusion junctions as those overlapping flanking + target exon
ei_events <- tumor_se_rmats %>%
  dplyr::filter(upstream_to_target_count >= 5) %>%
  dplyr::select(up_incl_jc) %>%
  dplyr::rename(junction = up_incl_jc) %>%
  bind_rows(tumor_se_rmats %>% 
              dplyr::filter(target_to_downstream_count >= 5) %>%
              dplyr::select(down_incl_jc) %>%
              dplyr::rename(junction = down_incl_jc)) %>%
  distinct() %>%
  dplyr::mutate(type = "EI") 

# define exon skipping events as those overlapping flanking junctions
es_events <- tumor_se_rmats %>%
  dplyr::filter(upstream_to_downstream_count >= 5) %>%
  dplyr::select(skip_jc) %>%
  dplyr::rename(junction = skip_jc) %>%
  distinct() %>%
  dplyr::mutate(type = "ES") 

# define intron retention events as those overlapping flanking exon and intron
ir_events <- tumor_ri_rmats %>%
  dplyr::filter(upstream_to_intron_count >= 5) %>%
  dplyr::select(up_incl_jc) %>%
  dplyr::rename(junction = up_incl_jc) %>%
  bind_rows(tumor_ri_rmats %>% 
              dplyr::filter(intron_to_downstream_count >= 5) %>%
              dplyr::select(down_incl_jc) %>%
              dplyr::rename(junction = down_incl_jc)) %>%
  distinct() %>%
  dplyr::mutate(type = "IR") 

# define intron exclusion events as those overlapping flanking exons 
ie_events <- tumor_ri_rmats %>%
  dplyr::filter(upstream_to_downstream_count >= 5) %>%
  dplyr::select(skip_jc) %>%
  dplyr::rename(junction = skip_jc) %>%
  distinct() %>%
  dplyr::mutate(type = "IE") 

# create merged junction annotation data frame
jc_annot <- ei_events %>%
  bind_rows(es_events, 
            ir_events,
            ie_events,
            a3_loss_events,
            a3_gain_events,
            a5_loss_events,
            a5_gain_events) %>%
  group_by(junction) %>%
  summarise(event_type = str_c(type, collapse = ";"))

# Merge tumor-specific junctions and targets 
merged_junction_df <- a5ss_junction_df %>%
  dplyr::select(geneSymbol, junction, boundary, strand, junction_count, n_low_ct_gtex) %>%
  bind_rows(a3ss_junction_df %>%
              dplyr::select(geneSymbol, junction, boundary, strand, 
                            junction_count, n_low_ct_gtex)) %>%
  bind_rows(se_junction_df %>%
              dplyr::select(geneSymbol, junction, boundary, strand, 
                            junction_count, n_low_ct_gtex)) %>%
  bind_rows(ri_junction_df %>%
              dplyr::select(geneSymbol, junction, boundary, strand, 
                            junction_count, n_low_ct_gtex)) %>%
  dplyr::rename(coordinates = junction,
                count = junction_count) %>%
  dplyr::mutate(type = "junction") %>%
  bind_rows(se_target_df %>%
              dplyr::select(geneSymbol, target, strand, target_count, 
                            n_low_ct_gtex) %>%
              dplyr::rename(coordinates = target,
                            count = target_count) %>%
              dplyr::mutate(type = "exon")) %>%
  bind_rows(ri_target_df %>%
              dplyr::select(geneSymbol, target, strand, target_count, 
                            n_low_ct_gtex) %>%
              dplyr::rename(coordinates = target,
                            count = target_count) %>%
              dplyr::mutate(type = "intron")) %>%
  distinct() %>%
  left_join(jc_annot,
            by = c("coordinates" = "junction")) %>%
  # define chr, start, end
  dplyr::mutate(chr = unlist(lapply(str_split(coordinates, ":"), function(x) x[1]))) %>%
  dplyr::mutate(start = case_when(
    type == "junction" ~ unlist(lapply(str_split(boundary, "_"), function(x) x[1])),
    type %in% c("exon", "intron") ~ unlist(lapply(str_split(coordinates, "[:|_]"), function(x) x[2]))
  ), end = case_when(
    type == "junction" ~ unlist(lapply(str_split(boundary, "_"), function(x) x[2])),
    type %in% c("exon", "intron") ~ unlist(lapply(str_split(coordinates, "[:|_]"), function(x) x[3]))
  ))

# determine if splice sites are annotated
merged_junction_df <- merged_junction_df %>%
  dplyr::mutate(start_annotated = NA_character_,
                end_annotated = NA_character_)
  
# load gencode GTF
gtf <- rtracklayer::import(gencode_gtf_file)
gtf_df <- data.frame(gtf)

# determine annotation of junction sites
for (i in seq_len(nrow(merged_junction_df))){
  
  type <- merged_junction_df$type[i]
  chr <- merged_junction_df$chr[i]
  start_pos <- ifelse(type == "junction",
                      as.numeric(merged_junction_df$start[i]),
                      as.numeric(merged_junction_df$start[i]) + 1)
  end_pos <- ifelse(type == "junction",
                as.numeric(merged_junction_df$end[i]) + 1,
                as.numeric(merged_junction_df$end[i]))
  
  start_annots <- gtf_df %>%
    dplyr::filter(seqnames == chr,
                  (start == start_pos | end == start_pos))
  
  end_annots <- gtf_df %>% 
    dplyr::filter(seqnames == chr,
                  (start == end_pos | end == end_pos))
  
  merged_junction_df$start_annotated[i] <- ifelse(nrow(start_annots) > 0,
                                               "Yes",
                                               "No")
  
  merged_junction_df$end_annotated[i] <- ifelse(nrow(end_annots) > 0,
                                               "Yes",
                                               "No")
  
}

# define whether junction contains any novel splice site
merged_junction_df <- merged_junction_df %>%
  dplyr::mutate(novel_splice_site = case_when(
    start_annotated == "No" | end_annotated == "No" ~ "Yes",
    TRUE ~ "No"
  ))

# save to output (CWL collects from working directory)
write_tsv(merged_junction_df,
          glue::glue("{sample_id}-tumor-specific-splice-products.tsv"))