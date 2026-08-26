"""
Overwrite a BEAN ReporterScreen object's per-guide "edits" layer with the
accessibility-adjusted scaling factor computed in Step 2, in place of BEAN's
own reporter-based edit calling.

BEAN's default edit-rate estimate for each guide is (edits layer) / (X_bcmatch
layer), derived from its own reporter allele calls. For this screen those
calls were replaced with the externally-computed per-guide scaling factor
(02_scale_counts_for_bean.R's edit_rates_one_per_guide.csv) so that BEAN's
model uses the same accessibility-adjusted predicted editing rate as the rest
of the pipeline: edits = scaling_factor * X_bcmatch. Guides missing from the
scaling-factor table fall back to using X_bcmatch directly (edit_rate = 1),
since there's no reporter-based signal to scale for them.

This runs between `bean qc` and `bean run` in 03_run_bean.sh.

Usage:
    python 03b_apply_manual_edit_rates.py <qc.h5ad> <edit_rates_one_per_guide.csv> <target_base_changes> <output.h5ad>
"""

import sys
import bean as be
import pandas as pd
import numpy as np

qc_h5ad_path = sys.argv[1]
edit_rates_path = sys.argv[2]
target_base_changes = sys.argv[3]  # e.g. "A>G"; required by the ReporterScreen object
out_path = sys.argv[4]

edit_rate_df = pd.read_csv(edit_rates_path, header=None, names=["grna_name", "edit_rate"])
bdata = be.read_h5ad(qc_h5ad_path)

bdata.uns["target_base_changes"] = target_base_changes

edit_rate_dict = dict(zip(edit_rate_df["grna_name"], edit_rate_df["edit_rate"]))
edit_rates = np.array([edit_rate_dict.get(g, np.nan) for g in bdata.obs.index])

missing_count = np.isnan(edit_rates).sum()
if missing_count > 0:
    print(f"Warning: {missing_count} guides have no scaling factor; using edit_rate = 1 (edits = X_bcmatch) for these")

X_bcmatch = bdata.layers["X_bcmatch"]
edit_rate_matrix = np.tile(edit_rates.reshape(-1, 1), (1, bdata.n_vars))
new_edits = np.where(np.isnan(edit_rate_matrix), X_bcmatch, edit_rate_matrix * X_bcmatch)

bdata.layers["edits"] = new_edits
bdata.obs["edit_rate"] = edit_rates

computed_edit_rate = bdata.layers["edits"] / (bdata.layers["X_bcmatch"] + 1e-10)
mean_edit_rates = np.nanmean(computed_edit_rate, axis=1)
print("Effective edit rate per guide, averaged across samples:")
print(f"  min={np.nanmin(mean_edit_rates):.4f} max={np.nanmax(mean_edit_rates):.4f} "
      f"mean={np.nanmean(mean_edit_rates):.4f} median={np.nanmedian(mean_edit_rates):.4f}")

bdata.write(out_path)
print(f"Wrote {out_path}")
