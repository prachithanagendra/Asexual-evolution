#!/usr/bin/env python3

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ==========================
# PATHS
# ==========================
INPUT_DIR = "/14_AD_AF_extracted"
BASE_OUTPUT = "/15_DPless30_AFdiff_heatmap"

DP_FILTER_DIR = os.path.join(BASE_OUTPUT, "01_DP30_filtered")
AF_FILTER_DIR = os.path.join(BASE_OUTPUT, "02_REF_AF_different")
HEATMAP_DIR = os.path.join(BASE_OUTPUT, "03_heatmaps")

os.makedirs(DP_FILTER_DIR, exist_ok=True)
os.makedirs(AF_FILTER_DIR, exist_ok=True)
os.makedirs(HEATMAP_DIR, exist_ok=True)

# ==========================
# PARAMETERS
# ==========================
ANCESTOR = "Ancestor2a"
DP_THRESHOLD = 30
AF_DIFF_THRESHOLD = 0.2

# ==========================
# PROCESS EACH SET
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

    # --------------------------
    # STEP 1: DP FILTER
    # Keep variant only if ALL samples have DP ≥ 30
    # --------------------------
    dp_columns = [f"{s}_DP" for s in samples]
    df_dp = df[df[dp_columns].min(axis=1) >= DP_THRESHOLD].copy()

    dp_out = os.path.join(DP_FILTER_DIR, tsv)
    df_dp.to_csv(dp_out, sep="\t", index=False)

    # --------------------------
    # STEP 2: REF_AF DIFFERENCE FILTER
    # Keep variant if at least one evolved sample differs
    # from ancestor by > 0.2
    # --------------------------
    ancestor_af = df_dp[f"{ANCESTOR}_REF_AF"]

    keep_mask = np.zeros(len(df_dp), dtype=bool)

    for s in samples:
        if s == ANCESTOR:
            continue
        diff = np.abs(df_dp[f"{s}_REF_AF"] - ancestor_af)
        keep_mask |= (diff > AF_DIFF_THRESHOLD)

    df_af = df_dp[keep_mask].copy()

    af_out = os.path.join(AF_FILTER_DIR, tsv)
    df_af.to_csv(af_out, sep="\t", index=False)

    # --------------------------
    # STEP 3: HEATMAP
    # --------------------------
    if len(df_af) == 0:
        print("  No variants left after filtering — skipping heatmap")
        continue

    af_matrix = df_af[[f"{s}_REF_AF" for s in samples]]

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
        HEATMAP_DIR,
        tsv.replace(".tsv", "_REF_AF_heatmap.png")
    )
    plt.savefig(heatmap_path, dpi=300)
    plt.close()

    print(f"  Variants after DP filter: {len(df_dp)}")
    print(f"  Variants after AF filter: {len(df_af)}")
    print("  Heatmap saved")

print("\nAll sets processed successfully.")

