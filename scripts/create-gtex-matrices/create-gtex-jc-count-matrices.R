## GTEx splice event junction & target count matrix generation 
##
## Ryan Corbett
##
## Nov 2025

# This script loads GTEx rMATS files, pull splice event junction and target coordinates, and calculates mean junction/target counts per GTEx subgroup

library(tidyverse)
library(qs2)
library(stringr)
library(data.table)

### Set directory paths
# Use current working directory (CWL stages files here)
data_dir <- getwd()
results_dir <- file.path("results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Set file paths (files are staged in current directory by CWL)
gtex_se_file <- "gtex-rmats_merged_raw_SE.qs2"
gtex_ri_file <- "gtex-rmats_merged_raw_RI.qs2"
gtex_a3ss_file <- "gtex-rmats_merged_raw_A3SS.qs2"
gtex_a5ss_file <- "gtex-rmats_merged_raw_A5SS.qs2"
hist_file <- "histologies.tsv"
gtex_read_ct_file <- file.path(data_dir, "gtex_input_read_counts.tsv")

# Load OPC hist and filter for GTEx RNA-seq
gtex_hist <- read_tsv(hist_file) %>%
  dplyr::filter(cohort == "GTEx",
                experimental_strategy == "RNA-Seq") %>%
  # only need BS ID, subgroup columns
  dplyr::select(Kids_First_Biospecimen_ID, gtex_subgroup)

# Load table of read counts used by rMATS per GTEx sample
gtex_read_cts <- read_tsv(gtex_read_ct_file) %>%
  dplyr::mutate(sample_id = unlist(lapply(strsplit(sample_id, "_"), function(x) x[[1]])))

# pull cell samples to be filtered out
gtex_cell_samples <- gtex_hist %>%
  dplyr::filter(grepl("Cells", gtex_subgroup)) %>%
  pull(Kids_First_Biospecimen_ID)

# Load rMATS SE results
se_df <- qs2::qs_read(gtex_se_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id)) %>%
  dplyr::filter(!sample_id %in% gtex_cell_samples)

# Define SE junction and target IDs, and select relevant columns 
se_df <- se_df %>%
  dplyr::mutate(
    up_incl_jc = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                       exonStart_0base, "-", exonEnd),
    down_incl_jc = str_c(chr, ":", exonStart_0base, "-", exonEnd, "_",
                         downstreamES, "-", downstreamEE),
    target = str_c(chr, ":", exonStart_0base, "_", exonEnd),
    skip_jc  = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                     downstreamES, "-", downstreamEE)
  ) %>%
  dplyr::select(sample_id, geneSymbol, up_incl_jc,
                down_incl_jc, target, skip_jc,
                upstream_to_target_count,
                target_to_downstream_count,
                target_count,
                upstream_to_downstream_count)

# convert to data table 
se_df <- as.data.table(se_df)

# get median target counts 
se_target_df <- se_df[
  , .(target_count = median(target_count, na.rm = TRUE)),
  by = .(sample_id, geneSymbol, target)
]

# convert read cts and hist to data tables
gtex_read_cts_dt <- as.data.table(gtex_read_cts)
gtex_hist_dt     <- as.data.table(gtex_hist)

# Merge histology info
se_target_df <- merge(
  se_target_df,
  gtex_hist_dt[, .(Kids_First_Biospecimen_ID, gtex_subgroup)],
  by.x = "sample_id",
  by.y = "Kids_First_Biospecimen_ID",
  all.x = TRUE
)

# Group and compute means
agg_se_target_df <- se_target_df[
  , .(mean_target_count = mean(target_count, na.rm = TRUE)),
  by = .(gtex_subgroup, target)
]

# Pivot wider
norm_se_target_mat <- dcast(
  agg_se_target_df,
  target ~ gtex_subgroup,
  value.var = "mean_target_count"
)

# Extract all junctions and pivot longer
se_junction_df <- rbindlist(list(
  se_df[, .(sample_id, geneSymbol,
            junction = up_incl_jc,
            junction_ct = upstream_to_target_count)],
  se_df[, .(sample_id, geneSymbol,
            junction = down_incl_jc,
            junction_ct = target_to_downstream_count)],
  se_df[, .(sample_id, geneSymbol,
            junction = skip_jc,
            junction_ct = upstream_to_downstream_count)]
))

# calculate median junction counts
se_junction_df <- se_junction_df[
  , .(junction_count = median(junction_ct)),
  by = .(sample_id, geneSymbol, junction)
]

# Merge histology info
se_junction_df <- merge(
  se_junction_df,
  gtex_hist_dt[, .(Kids_First_Biospecimen_ID, gtex_subgroup)],
  by.x = "sample_id",
  by.y = "Kids_First_Biospecimen_ID",
  all.x = TRUE
)

# Group and compute means
agg_se_junction_df <- se_junction_df[
  , .(mean_junction_count = mean(junction_count, na.rm = TRUE)),
  by = .(gtex_subgroup, junction)
]

# Wide format
norm_se_junction_mat <- dcast(
  agg_se_junction_df,
  junction ~ gtex_subgroup,
  value.var = "mean_junction_count"
)

# save PSI, target, and junction count matrices
qs2::qs_save(norm_se_target_mat,
             file.path(results_dir, "gtex-se-target-mat-raw.qs2"))
qs2::qs_save(norm_se_junction_mat,
             file.path(results_dir, "gtex-se-junction-mat-raw.qs2"))

rm(se_df, se_junction_df, se_target_df, norm_se_target_mat, norm_se_junction_mat)

# Load retained intron (RI) rMATS results, update sample ID, and filter out cell samples
ri_df <- qs2::qs_read(gtex_ri_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id)) %>%
  dplyr::filter(!sample_id %in% gtex_cell_samples)

# define columns specifying junction coordinates and retain only relevant columns
ri_df <- ri_df %>%
  dplyr::mutate(
    up_incl_jc = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                       upstreamEE, "-", downstreamES),
    
    down_incl_jc = str_c(chr, ":", upstreamEE, "-", downstreamES, "_",
                         downstreamES, "-", downstreamEE),
    
    target = str_c(chr, ":", upstreamEE, "_", downstreamES),
    
    skip_jc = str_c(chr, ":", upstreamES, "-", upstreamEE, "_",
                    downstreamES, "-", downstreamEE)
  ) %>%
  dplyr::rename(target_count = intron_count) %>%
  dplyr::select(sample_id, geneSymbol, up_incl_jc,
                down_incl_jc, target, skip_jc,
                upstream_to_intron_count,
                intron_to_downstream_count,
                target_count,
                upstream_to_downstream_count)

# convert to data table 
ri_df <- as.data.table(ri_df)

# get median intron counts 
ri_target_df <- ri_df[
  , .(target_count = median(target_count)),
  by = .(sample_id, geneSymbol, target)
]

# Merge histology info
ri_target_df <- merge(
  ri_target_df,
  gtex_hist_dt[, .(Kids_First_Biospecimen_ID, gtex_subgroup)],
  by.x = "sample_id",
  by.y = "Kids_First_Biospecimen_ID",
  all.x = TRUE
)

# Group and compute means
agg_ri_target_df <- ri_target_df[
  , .(mean_target_count = mean(target_count, na.rm = TRUE)),
  by = .(gtex_subgroup, target)
]

# Wide format
norm_ri_target_mat <- dcast(
  agg_ri_target_df,
  target ~ gtex_subgroup,
  value.var = "mean_target_count"
)

# Build the long-form junction table, merging data from:
# upstream-intron 
# intron-downstream
# upstream-downstream junctions
ri_junction_df <- rbindlist(list(
  ri_df[, .(sample_id, geneSymbol,
            junction = up_incl_jc,
            junction_ct = upstream_to_intron_count)],
  ri_df[, .(sample_id, geneSymbol,
            junction = down_incl_jc,
            junction_ct = intron_to_downstream_count)],
  ri_df[, .(sample_id, geneSymbol,
            junction = skip_jc,
            junction_ct = upstream_to_downstream_count)]
))

# calculate median junction counts
ri_junction_df <- ri_junction_df[
  , .(junction_count = median(junction_ct)),
  by = .(sample_id, geneSymbol, junction)
]

# Merge histology info
ri_junction_df <- merge(
  ri_junction_df,
  gtex_hist_dt[, .(Kids_First_Biospecimen_ID, gtex_subgroup)],
  by.x = "sample_id",
  by.y = "Kids_First_Biospecimen_ID",
  all.x = TRUE
)

# Group and compute means
agg_ri_junction_df <- ri_junction_df[
  , .(mean_junction_count = mean(junction_count, na.rm = TRUE)),
  by = .(gtex_subgroup, junction)
]

# Wide format
norm_ri_junction_mat <- dcast(
  agg_ri_junction_df,
  junction ~ gtex_subgroup,
  value.var = "mean_junction_count"
)

# save PSI, target, and junction count matrices
qs2::qs_save(norm_ri_target_mat,
             file.path(results_dir, "gtex-ri-target-mat-raw.qs2"))
qs2::qs_save(norm_ri_junction_mat,
             file.path(results_dir, "gtex-ri-junction-mat-raw.qs2"))

rm(ri_df, ri_junction_df, ri_target_df, norm_ri_target_mat, norm_ri_junction_mat)

# A3SS events, modify sample ID and filter out cell samples 
a3ss_df <- qs2::qs_read(gtex_a3ss_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id)) %>%
  dplyr::filter(!sample_id %in% gtex_cell_samples)

# define junction coordinates and filter for relevant columns
a3ss_df <- a3ss_df %>%
  dplyr::mutate(
    long_incl_jc = str_c(chr, ":",
                         flankingES, "-", flankingEE, "_",
                         longExonStart_0base, "-", longExonEnd),
    
    short_incl_jc = str_c(chr, ":",
                          flankingES, "-", flankingEE, "_",
                          shortES, "-", shortEE)
  ) %>%
  dplyr::select(sample_id, geneSymbol, long_incl_jc,
                short_incl_jc,
                long_to_flanking_count,
                short_to_flanking_count)

# convert to data table
a3ss_df <- as.data.table(a3ss_df)

# pivot longer for single row per unique sample & junction 
a3ss_junction_df <- rbindlist(list(
  a3ss_df[, .(sample_id, geneSymbol,
              junction = long_incl_jc,
              junction_ct = long_to_flanking_count)],
  a3ss_df[, .(sample_id, geneSymbol,
              junction = short_incl_jc,
              junction_ct = short_to_flanking_count)]
))

# calculate median junction counts when duplicated 
a3ss_junction_df <- a3ss_junction_df[
  , .(junction_count = median(junction_ct)),
  by = .(sample_id, geneSymbol, junction)
]

# Merge histology info
a3ss_junction_df <- merge(
  a3ss_junction_df,
  gtex_hist_dt[, .(Kids_First_Biospecimen_ID, gtex_subgroup)],
  by.x = "sample_id",
  by.y = "Kids_First_Biospecimen_ID",
  all.x = TRUE
)

# Group and compute means
agg_a3ss_junction_df <- a3ss_junction_df[
  , .(mean_junction_count = mean(junction_count, na.rm = TRUE)),
  by = .(gtex_subgroup, junction)
]

# Wide format
norm_a3ss_junction_mat <- dcast(
  agg_a3ss_junction_df,
  junction ~ gtex_subgroup,
  value.var = "mean_junction_count"
)

# save junction count matrix
qs2::qs_save(norm_a3ss_junction_mat,
             file.path(results_dir, "gtex-a3ss-junction-mat-raw.qs2"))

rm(a3ss_df, a3ss_junction_df, norm_a3ss_junction_mat)

# A5SS events
a5ss_df <- qs2::qs_read(gtex_a5ss_file) %>%
  dplyr::mutate(sample_id = sub("_.*", "", sample_id)) %>%
  dplyr::filter(!sample_id %in% gtex_cell_samples)

# define junction coordinates and filter for relevant columns 
a5ss_df <- a5ss_df %>%
  dplyr::mutate(
    long_incl_jc = str_c(chr, ":",
                         longExonStart_0base, "-", longExonEnd, "_",
                         flankingES, "-", flankingEE),
    
    short_incl_jc = str_c(chr, ":",
                          shortES, "-", shortEE, "_",
                          flankingES, "-", flankingEE)
  ) %>%
  dplyr::select(sample_id, geneSymbol, long_incl_jc,
                short_incl_jc,
                long_to_flanking_count,
                short_to_flanking_count)

# convert to data.table
a5ss_df <- as.data.table(a5ss_df)

# pivot longer for single row per unique sample & junction 
a5ss_junction_df <- rbindlist(list(
  a5ss_df[, .(sample_id, geneSymbol,
              junction = long_incl_jc,
              junction_ct = long_to_flanking_count)],
  a5ss_df[, .(sample_id, geneSymbol,
              junction = short_incl_jc,
              junction_ct = short_to_flanking_count)]
))

# calculate median junction counts
a5ss_junction_df <- a5ss_junction_df[
  , .(junction_count = median(junction_ct)),
  by = .(sample_id, geneSymbol, junction)
]

# Merge histology info
a5ss_junction_df <- merge(
  a5ss_junction_df,
  gtex_hist_dt[, .(Kids_First_Biospecimen_ID, gtex_subgroup)],
  by.x = "sample_id",
  by.y = "Kids_First_Biospecimen_ID",
  all.x = TRUE
)

# Group and compute means
agg_a5ss_junction_df <- a5ss_junction_df[
  , .(mean_junction_count = mean(junction_count, na.rm = TRUE)),
  by = .(gtex_subgroup, junction)
]

# Wide format
norm_a5ss_junction_mat <- dcast(
  agg_a5ss_junction_df,
  junction ~ gtex_subgroup,
  value.var = "mean_junction_count"
)

# save junction count matrix
qs2::qs_save(norm_a5ss_junction_mat,
             file.path(results_dir, "gtex-a5ss-junction-mat-raw.qs2"))

rm(a5ss_df, a5ss_junction_df, norm_a5ss_junction_mat)
