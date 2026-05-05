#!/usr/bin/env python3

import pandas as pd

# ==========================================================
# USER INPUTS
# ==========================================================

# Main TSV file (the one that will be filtered)
INPUT_FILE1 = "/15_BiallelicSNP_analysis/2_DP20_AF_diff_heatmap/02_REF_AF_different/B12_B24.tsv"

# Fixed variants TSV file
INPUT_FILE2 = "/14_SNP_analysis/5_fixed_annotated_merged/gal6_B_annotation_AF_merged.tsv"

# Output TSV file
OUTPUT_FILE = "/15_BiallelicSNP_analysis/3_wo_fixed_variants/gal6_B_present.tsv"

# ==========================================================
# LOAD FILES
# ==========================================================

df1 = pd.read_csv(INPUT_FILE1, sep="\t")
df2 = pd.read_csv(INPUT_FILE2, sep="\t")

# ==========================================================
# CREATE VARIANT KEY (CHROM, POS, REF, ALT)
# ==========================================================

variants_file2 = set(
    zip(df2["CHROM"], df2["POS"], df2["REF"], df2["ALT"])
)

# ==========================================================
# FILTER FILE1
# ==========================================================

mask = ~df1.apply(
    lambda row: (row["CHROM"], row["POS"], row["REF"], row["ALT"]) in variants_file2,
    axis=1
)

filtered_df = df1[mask]

# ==========================================================
# SAVE OUTPUT
# ==========================================================

filtered_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print("Filtering complete.")
print(f"Original rows in file1: {len(df1)}")
print(f"Rows after removing fixed variants: {len(filtered_df)}")
print(f"Output written to: {OUTPUT_FILE}")
