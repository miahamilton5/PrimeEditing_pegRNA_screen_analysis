#!/usr/bin/env Rscript
# Plot the distribution of per-pegRNA reporter editing rates (average percent
# editing across bulk replicates, from 01_generate_count_table.R), split into
# within-screen quartiles.
#
# Usage:
#   Rscript plot_reporter_editing_density.R <average_editing.csv> <output_plot.png>

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggridges))

args <- commandArgs(trailingOnly = TRUE)
editing_path <- args[1]
out_path <- args[2]

editing <- fread(editing_path) %>%
  mutate(quartile = ntile(average_percent_editing, 4))

cat("median editing:", median(editing$average_percent_editing, na.rm = TRUE), "\n")
cat("mean editing:", mean(editing$average_percent_editing, na.rm = TRUE), "\n")

p <- ggplot(editing, aes(x = average_percent_editing, y = 0, fill = factor(quartile))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE, quantiles = 4, quantile_lines = TRUE) +
  scale_x_continuous(limits = c(-15, 100)) +
  theme_classic() +
  labs(x = "Reporter editing %", y = NULL, fill = "Quartile")

ggsave(out_path, p, width = 6, height = 3, dpi = 200)
