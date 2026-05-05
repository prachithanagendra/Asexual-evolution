#!/usr/bin/env python3

import os
import pandas as pd

# === Input folders ===
folder_del = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/deletion_annotated/unique_genes_del"
folder_dup = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/duplication_annotated/unique_genes_dup"

# === Output files ===
output_del = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/RD_deletion_gene_list.tsv"
output_dup = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/RD_duplication_gene_list.tsv"


def compile_gene_lists(folder, suffix=""):

    sample_dict = {}
    max_len = 0
    missing_samples = []

    for i in range(1, 49):

        sample_name = f"A{i}"
        file_path = os.path.join(folder, f"{sample_name}_unique_gene_names{suffix}.txt")

        if os.path.exists(file_path):
            with open(file_path) as f:
                genes = [line.strip() for line in f if line.strip()]
        else:
            genes = []   # empty column if file missing
            missing_samples.append(sample_name)

        sample_dict[sample_name] = genes
        max_len = max(max_len, len(genes))

    # Pad lists so all columns have equal length
    for sample in sample_dict:
        sample_dict[sample] += [""] * (max_len - len(sample_dict[sample]))

    df = pd.DataFrame(sample_dict)

    if missing_samples:
        print(f"Missing files for: {', '.join(missing_samples)}")

    return df


# === Compile deletion and duplication gene tables ===
df_del = compile_gene_lists(folder_del, suffix="_del")
df_dup = compile_gene_lists(folder_dup, suffix="_dup")

# === Save files ===
df_del.to_csv(output_del, sep="\t", index=False)
df_dup.to_csv(output_dup, sep="\t", index=False)

print(f"Deletion gene list saved: {output_del}")
print(f"Duplication gene list saved: {output_dup}")
