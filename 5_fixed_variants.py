#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# PATHS
# ==========================
INPUT_DIR = "/15_DPless30_AFdiff_heatmap/02_REF_AF_different"

OUTPUT_BASE = "/16_fixed_variants"
OUTPUT_TSV_DIR = os.path.join(OUTPUT_BASE, "tsv")
OUTPUT_HEATMAP_DIR = os.path.join(OUTPUT_BASE, "heatmaps")

os.makedirs(OUTPUT_TSV_DIR, exist_ok=True)
os.makedirs(OUTPUT_HEATMAP_DIR, exist_ok=True)

# ==========================
# PARAMETERS
# ==========================
ANCESTOR = "Ancestor2a"

LOW = (0.0, 0.2)
MID = (0.4, 0.6)
HIGH = (0.8, 1.0)

# ==========================
# Helper function
# ==========================
def in_range(x, r):
    return (x >= r[0]) & (x <= r[1])

# ==========================
# PROCESS EACH FILE
# ==========================
tsv_files = sorted([f for f in os.listdir(INPUT_DIR) if f.endswith(".tsv")])

for tsv in tsv_files:
    print(f"\nProcessing {tsv}")
    df = pd.read_csv(os.path.join(INPUT_DIR, tsv), sep="\t")

    # --------------------------
    # Identify samples (ordered)
    # --------------------------
    dp_cols = [c for c in df.columns if c.endswith("_DP")]
    samples = [c.replace("_DP", "") for c in dp_cols]

    # Ancestor first, then numeric order
    samples = sorted(
        samples,
        key=lambda x: (x != ANCESTOR, int(x[1:]) if x != ANCESTOR else -1)
    )

    ancestor_af = df[f"{ANCESTOR}_REF_AF"]

    # --------------------------
    # Categorical REF_AF filter
    # --------------------------
    keep_mask = np.zeros(len(df), dtype=bool)

    for s in samples:
        if s == ANCESTOR:
            continue

        sample_af = df[f"{s}_REF_AF"]

        # Rule 1: Anc low -> sample mid/high
        rule1 = (
            in_range(ancestor_af, LOW) &
            (in_range(sample_af, MID) | in_range(sample_af, HIGH))
        )

        # Rule 2: Anc high -> sample low/mid
        rule2 = (
            in_range(ancestor_af, HIGH) &
            (in_range(sample_af, LOW) | in_range(sample_af, MID))
        )

        # Rule 3: Anc mid -> sample low/high
        rule3 = (
            in_range(ancestor_af, MID) &
            (in_range(sample_af, LOW) | in_range(sample_af, HIGH))
        )

        keep_mask |= (rule1 | rule2 | rule3)

    df_out = df[keep_mask].copy()

    # --------------------------
    # Write filtered TSV
    # --------------------------
    out_tsv = os.path.join(OUTPUT_TSV_DIR, tsv)
    df_out.to_csv(out_tsv, sep="\t", index=False)

    print(f"  Variants kept: {len(df_out)}")

    # --------------------------
    # Heatmap (REF_AF)
    # --------------------------
    if len(df_out) == 0:
        print("  No variants left — skipping heatmap")
        continue

    af_matrix = df_out[[f"{s}_REF_AF" for s in samples]]

    plt.figure(figsize=(1.2 * len(samples), 0.15 * len(af_matrix) + 3))
    plt.imshow(af_matrix, aspect="auto", interpolation="nearest")
    plt.colorbar(label="REF_AF")

    plt.xticks(
        ticks=range(len(samples)),
        labels=samples,
        rotation=45,
        ha="right"
    )
    plt.yticks([])

    plt.title(tsv.replace(".tsv", ""))
    plt.tight_layout()

    heatmap_path = os.path.join(
        OUTPUT_HEATMAP_DIR,
        tsv.replace(".tsv", "_REF_AF_category_heatmap.png")
    )

    plt.savefig(heatmap_path, dpi=300)
    plt.close()

    print("  Heatmap saved")

print("\nAll categorical filtering and heatmaps complete.")

