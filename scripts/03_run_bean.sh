#!/usr/bin/env bash
# Build a BEAN ReporterScreen object from the count table, run BEAN's QC step,
# and quantify variant effect sizes.
#
# BEAN (Ryu et al. 2024, Nat Genet; https://github.com/pinellolab/crispr-bean)
# is run in "variant" library-design mode, since several pegRNAs tile each
# targeted variant and bystander edits outside the intended edit are not of
# interest here (unlike a coding-sequence tiling screen).
#
# The accessibility-adjusted editing-rate scaling from Step 2 is incorporated
# into the ReporterScreen object's per-guide metadata before this run (via a
# manual correction step after `bean qc`, not included in this repo) rather
# than passed as `bean run` flags -- the actual `bean run` call for this
# screen used no extra options beyond `-o`.
#
# Usage:
#   ./03_run_bean.sh <pegRNA_library.csv> <sample_info.csv> <edited_counts.csv> <edits.csv> <bcmatch.csv> <output_dir>

set -euo pipefail

LIBRARY_CSV="$1"      # gRNA/pegRNA info table: name, target, target_group, scaling, ...
SAMPLE_INFO_CSV="$2"  # sample info table: condition (bulk/b20/t20), replicate
COUNTS_CSV="$3"       # oligo x sample raw count matrix
EDITS_CSV="$4"        # oligo x sample edited-pseudo-count matrix from 02_scale_counts_for_bean.R
BCMATCH_CSV="$5"      # oligo x sample barcode-matched count matrix
OUT_DIR="$6"

mkdir -p "$OUT_DIR"

bean create-screen \
  -e "$EDITS_CSV" \
  -b "$BCMATCH_CSV" \
  "$LIBRARY_CSV" \
  "$SAMPLE_INFO_CSV" \
  "$COUNTS_CSV" \
  -o "$OUT_DIR/bean_counts"

bean qc \
  "$OUT_DIR/bean_counts.h5ad" \
  -o "$OUT_DIR/bean_counts_QC.h5ad" \
  -r "$OUT_DIR/qc_report" \
  --dont-recalculate-edits

# At this point, per-guide editing rates and target-base edit calls in
# bean_counts_QC.h5ad are manually reviewed and corrected as needed (this
# step is not scripted here); the corrected object is what's passed to
# `bean run` below.

bean run sorting variant \
  "$OUT_DIR/bean_counts_QC.h5ad" \
  -o "$OUT_DIR"

echo "BEAN element-level results: $OUT_DIR/bean_element_result.MixtureNormal.csv"
