#!/usr/bin/env python3

import os
import pandas as pd

# ==========================================================
# USER INPUT
# ==========================================================

BASE_INPUT_DIR  = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/2_cnv_calls"
BASE_OUTPUT_DIR = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/3_merged_calls"

os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

# ==========================================================
# SAMPLE GROUP DEFINITIONS
# ==========================================================

groups = {
    "glu_300":  range(1, 7),
    "gal_300":  range(7, 13),
    "glu_600":  range(13, 19),
    "gal_600":  range(19, 25),
    "glu_900":  range(25, 31),
    "gal_900":  range(31, 37),
    "glu_1200": range(37, 43),
    "gal_1200": range(43, 49),
}

# ==========================================================
# LOOP OVER GROUPS
# ==========================================================

for group_name, sample_range in groups.items():

    print(f"\nProcessing {group_name}...")

    merged_df = None

    for i in sample_range:
        sample = f"A{i}"
        filepath = os.path.join(BASE_INPUT_DIR, f"{sample}_sorted1_cnv_calls.tsv")

        if not os.path.exists(filepath):
            print(f"⚠ Missing file: {filepath}")
            continue

        df = pd.read_csv(
            filepath,
            sep="\t",
            usecols=["chrom", "window", "end", "CNV_state"]
        )

        df = df.rename(columns={"CNV_state": sample})

        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(
                merged_df,
                df,
                on=["chrom", "window", "end"],
                how="outer"
            )

    if merged_df is not None:
        merged_df = merged_df.sort_values(by=["chrom", "window"])

        output_path = os.path.join(
            BASE_OUTPUT_DIR,
            f"merged_CNV_state_{group_name}.tsv"
        )

        merged_df.to_csv(output_path, sep="\t", index=False)

        print(f" Saved: {output_path}")
    else:
        print(f" No files merged for {group_name}")

print("\n🎉 All groups processed.")

