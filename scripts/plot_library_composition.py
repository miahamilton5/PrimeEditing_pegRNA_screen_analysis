"""
Plot the variant-class composition of a pegRNA library: number of pegRNAs
and number of unique targeted variants per class (SNV, synthetic indel,
TSS synthetic indel, common indel, matched control, non-targeting).

Usage:
    python plot_library_composition.py <pegRNA_library.csv> <output_plot.png>
"""

import sys
import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CLASS_ORDER = ["snv", "synthetic_indel", "TSS_synthetic_indel", "common_insertion", "matched_control", "non-targeting"]
CLASS_COLORS = {
    "snv": "#6699CC",
    "synthetic_indel": "#66CC99",
    "TSS_synthetic_indel": "#EC5F67",
    "common_insertion": "#CC99CC",
    "matched_control": "#999999",
    "non-targeting": "#CCCCCC",
}


def main():
    lib_path = sys.argv[1]
    out_path = sys.argv[2]

    pegrnas_per_class = defaultdict(int)
    variants_per_class = defaultdict(set)
    with open(lib_path) as f:
        for row in csv.DictReader(f):
            cls = row["class"]
            pegrnas_per_class[cls] += 1
            variants_per_class[cls].add(row["target"])

    classes = [c for c in CLASS_ORDER if c in pegrnas_per_class]
    pegrna_counts = [pegrnas_per_class[c] for c in classes]
    variant_counts = [len(variants_per_class[c]) for c in classes]
    colors = [CLASS_COLORS[c] for c in classes]

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
    axes[0].bar(classes, pegrna_counts, color=colors, edgecolor="white")
    axes[0].set_ylabel("Number of pegRNAs")
    axes[0].set_title("pegRNAs per class")
    axes[0].tick_params(axis="x", rotation=45)

    axes[1].bar(classes, variant_counts, color=colors, edgecolor="white")
    axes[1].set_ylabel("Number of unique targets")
    axes[1].set_title("Targets per class")
    axes[1].tick_params(axis="x", rotation=45)

    for ax in axes:
        for label in ax.get_xticklabels():
            label.set_ha("right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)

    total_pegrnas = sum(pegrna_counts)
    total_variants = sum(variant_counts)
    print(f"total pegRNAs: {total_pegrnas}, total unique targets: {total_variants}")
    for c, p, v in zip(classes, pegrna_counts, variant_counts):
        print(f"  {c}: {p} pegRNAs, {v} targets")


if __name__ == "__main__":
    main()
