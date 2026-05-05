#!/usr/bin/env python3

import os
import pandas as pd

folder_del = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/deletion_annotated/unique_genes_del"
folder_dup = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/duplication_annotated/unique_genes_dup"

total_genes = 6281

samples = []
count_del = []
count_dup = []

missing_del = []
missing_dup = []

for i in range(1, 49):

    sample = f"A{i}"

    del_file = os.path.join(folder_del, f"{sample}_unique_gene_names_del.txt")
    dup_file = os.path.join(folder_dup, f"{sample}_unique_gene_names_dup.txt")

    # deletion count
    if os.path.exists(del_file):
        with open(del_file) as f:
            del_count = sum(1 for _ in f)
    else:
        del_count = 0
        missing_del.append(sample)

    # duplication count
    if os.path.exists(dup_file):
        with open(dup_file) as f:
            dup_count = sum(1 for _ in f)
    else:
        dup_count = 0
        missing_dup.append(sample)

    samples.append(sample)
    count_del.append(del_count)
    count_dup.append(dup_count)

# Create dataframe
df = pd.DataFrame({
    "sample": samples,
    "count_del": count_del,
    "count_dup": count_dup
})

# Calculate percentages
df["del_ratio"] = (df["count_del"] / total_genes) * 100
df["dup_ratio"] = (df["count_dup"] / total_genes) * 100

# Save TSV
output_file = "/home1/NGF_SK1_reference/NGF_SK1_ref_analysis/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/CNV_gene_counts_summary.tsv"
df.to_csv(output_file, sep="\t", index=False)

print(f"Summary file saved as {output_file}")

# Only print missing files to terminal
if missing_del:
    print("\nDeletion file missing for:", ", ".join(missing_del))

if missing_dup:
    print("\nDuplication file missing for:", ", ".join(missing_dup))
