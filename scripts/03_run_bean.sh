#!/usr/bin/env bash
# Build a BEAN ReporterScreen object from the count table and run BEAN's
# variant-effect model.
#
# BEAN (Ryu et al. 2024, Nat Genet; https://github.com/pinellolab/crispr-bean)
# is run in "variant" library-design mode, since several pegRNAs tile each
# targeted variant and bystander edits outside the intended edit are not of
# interest here (unlike a coding-sequence tiling screen). Per the Methods:
# "BEAN was run with default parameters, pegRNAs grouped by edit, and scaled
# by ATAC score (minimum scaling of 0.1). Matched control pegRNAs were
# scaled equally to the corresponding targeting pegRNA."
#
# Usage:
#   ./03_run_bean.sh <pegRNA_library_with_scaling.csv> <sample_info.csv> <edited_counts.csv> <output_dir>

set -euo pipefail

LIBRARY_CSV="$1"      # gRNA/pegRNA info table: name, target, target_group, scaling, ...
SAMPLE_INFO_CSV="$2"  # sample info table: condition (bulk/b20/t20), replicate
EDITED_COUNTS_CSV="$3" # oligo x sample edited-pseudo-count matrix from 02_scale_counts_for_bean.R
OUT_DIR="$4"

mkdir -p "$OUT_DIR"

bean create-screen \
  "$LIBRARY_CSV" \
  "$SAMPLE_INFO_CSV" \
  "$EDITED_COUNTS_CSV" \
  -o "$OUT_DIR/bean_counts"

bean run sorting variant \
  "$OUT_DIR/bean_counts.h5ad" \
  --scale-by-acc \
  --acc-col scaling \
  --guide-activity-col average_percent_editing \
  --fit-negctrl \
  -o "$OUT_DIR"

echo "BEAN element-level results: $OUT_DIR/bean_element_result.MixtureNormal.csv"
