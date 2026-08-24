#!/usr/bin/env Rscript
# Compute allelic (targeting vs. matched-control) effect sizes and z-scores
# from a BEAN MixtureNormal element-level result table, and plot a volcano.
#
# BEAN reports one posterior effect size (mu) and standard deviation (mu_sd)
# per pegRNA-defined target, including both targeting pegRNAs and their
# matched controls (which install the reference allele at the same site to
# control for steric effects of the pegRNA/reporter construct itself, see
# docs/manuscript_figures/pegRNA_matched_control_schematic.png). To isolate
# the allele-specific effect of each variant, we subtract each targeting
# pegRNA's matched-control effect size from its own, and propagate the two
# posterior standard deviations into a z-score:
#
#   allelic_mu      = mu_targeting - mu_matched
#   allelic_z_score = (mu_targeting - mu_matched) / sqrt(mu_sd_targeting^2 + mu_sd_matched^2)
#
# Usage:
#   Rscript 04_analyze_bean_results.R <bean_element_result.MixtureNormal.csv> <grouped_matched_key.csv> <output_table.csv>

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
bean_result_path <- args[1]
matched_key_path <- args[2]
out_path <- args[3]

Z_THRESHOLD <- 1.98  # |Z-score| threshold used to call significant hits

bean <- fread(bean_result_path)
matched_key <- fread(matched_key_path)

targeting <- bean %>%
  filter(target_group != "matched_control") %>%
  select(target, n_guides, edit_rate_mean, mu, mu_sd, mu_z) %>%
  inner_join(matched_key, by = c("target" = "id"))

matched <- bean %>%
  filter(target_group == "matched_control") %>%
  select(target, mu, mu_sd) %>%
  rename(mu_matched = mu, mu_sd_matched = mu_sd)

results <- targeting %>%
  inner_join(matched, by = c("matched_control_variant" = "target")) %>%
  rename(mu_targeting = mu, mu_sd_targeting = mu_sd) %>%
  mutate(
    allelic_mu = mu_targeting - mu_matched,
    allelic_z_score = (mu_targeting - mu_matched) / sqrt(mu_sd_targeting^2 + mu_sd_matched^2),
    variant_type = case_when(
      grepl("TSS_synthetic_indel", target) ~ "5 bp TSS deletion",
      grepl("synthetic_indel", target) ~ "5 bp synthetic deletion",
      grepl("snv", target) ~ "SNP",
      TRUE ~ "other"
    ),
    significant = abs(allelic_z_score) >= Z_THRESHOLD
  )

fwrite(results, out_path)

cat("Variants tested:", nrow(results), "\n")
cat("Significant hits (|Z| >=", Z_THRESHOLD, "):", sum(results$significant), "\n")
print(results %>% filter(significant) %>% count(variant_type))

ggplot(results, aes(x = allelic_mu, y = abs(allelic_z_score), color = variant_type)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = Z_THRESHOLD, linetype = "dashed") +
  scale_color_manual(values = c(
    "SNP" = "#9B59B6",
    "5 bp synthetic deletion" = "#6699CC",
    "5 bp TSS deletion" = "#EC5F67",
    "other" = "lightgray"
  )) +
  theme_classic() +
  labs(x = "Allelic effect size (mu, targeting - matched control)", y = "|Z-score|", color = "Variant type")
