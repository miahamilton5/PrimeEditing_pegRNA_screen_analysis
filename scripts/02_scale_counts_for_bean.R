#!/usr/bin/env Rscript
# Scale a raw pegRNA count table into the "edited pseudo-count" and
# per-guide editing-rate inputs BEAN uses to estimate variant effect sizes.
#
# BEAN needs to know, for each pegRNA, what fraction of sequencing reads in
# a given bin actually carry the intended edit (reporter-measured editing
# efficiency is only a proxy for this at the endogenous locus, and is itself
# noisy for lowly-accessible loci). We correct the raw reporter editing rate
# for local chromatin accessibility using a linear model fit against ATAC-seq
# signal (see Methods; the same relationship used to predict endogenous
# editing from reporter editing in Figure 4J), then convert that into a
# per-guide scaling factor with a floor so that guides with very low
# predicted editing aren't scaled to near-zero weight.
#
# Usage:
#   Rscript 02_scale_counts_for_bean.R <counts.csv> <average_editing.csv> <atac_signal.csv> <pegRNA_library.csv> <output_dir>
#
# atac_signal.csv: columns oligo_id, avg_score_5bp_rollavg (mean ATAC-seq
# signal in a 50-bp window around each pegRNA's target site)

suppressMessages(library(dplyr))
suppressMessages(library(data.table))

MIN_SCALING <- 0.1

args <- commandArgs(trailingOnly = TRUE)
counts_path <- args[1]
editing_path <- args[2]
atac_path <- args[3]
library_path <- args[4]
out_dir <- args[5]

counts <- fread(counts_path)
editing <- fread(editing_path)
atac <- fread(atac_path)
lib <- fread(library_path)

sample_cols <- setdiff(colnames(counts), "oligo_id")

scaled <- editing %>%
  left_join(atac, by = "oligo_id") %>%
  left_join(lib %>% select(name, class), by = c("oligo_id" = "name")) %>%
  mutate(
    # ATAC-adjusted predicted endogenous editing rate, clamped to [0, 100]
    predicted_endo_editing = average_percent_editing * ((log10(avg_score_5bp_rollavg) * 1.5) - 0.34),
    predicted_endo_editing = pmin(pmax(predicted_endo_editing, 0), 100),
    scaling = case_when(
      predicted_endo_editing <= 10 ~ MIN_SCALING,
      class %in% c("matched_control", "non-targeting") ~ 1,
      TRUE ~ predicted_endo_editing / 100
    )
  )

# edited pseudo-count table: multiply each bin/rep raw count by the
# oligo's reporter-measured editing fraction, so BEAN's input read counts
# approximate the number of *edited* molecules observed in each bin
edited_counts <- counts %>%
  left_join(editing %>% select(oligo_id, average_percent_editing), by = "oligo_id") %>%
  mutate(across(all_of(sample_cols), ~ .x * average_percent_editing / 100)) %>%
  select(-average_percent_editing)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(edited_counts %>% arrange(oligo_id), file.path(out_dir, "edited_counts.csv"))

# BEAN's --guide-info per-guide scaling-factor file: one scaling value per
# oligo, repeated across all bin/rep columns to match the count matrix shape
scaling_matrix <- scaled %>%
  arrange(oligo_id) %>%
  select(oligo_id, scaling) %>%
  { bind_cols(select(., oligo_id), as.data.frame(matrix(rep(.$scaling, length(sample_cols)), ncol = length(sample_cols)))) }
fwrite(scaling_matrix, file.path(out_dir, "edit_rates.csv"), col.names = FALSE)
fwrite(scaled %>% select(oligo_id, scaling) %>% arrange(oligo_id), file.path(out_dir, "edit_rates_one_per_guide.csv"), col.names = FALSE)

cat("pegRNAs scaled:", nrow(scaled), "\n")
cat("pegRNAs at the", MIN_SCALING, "scaling floor:", sum(scaled$scaling == MIN_SCALING), "\n")
