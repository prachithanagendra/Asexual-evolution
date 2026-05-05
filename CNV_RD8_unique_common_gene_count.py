#!/usr/bin/env python3
import pandas as pd

# === Input/Output paths ===
input_file = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/RD_deletion_gene_list.tsv"   # change this to your file path
output_file = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/RD_deletion_unique_gene_counts.tsv"

# === Read file ===
df = pd.read_csv(input_file, sep="\t")

# Define the 12 sets with labels (time order: 300,600,900,1200)
sets = {
    "glu1": ["A1","A13","A25","A37"],
    "glu2": ["A2","A14","A26","A38"],
    "glu3": ["A3","A15","A27","A39"],
    "glu4": ["A4","A16","A28","A40"],
    "glu5": ["A5","A17","A29","A41"],
    "glu6": ["A6","A18","A30","A42"],
    "gal1": ["A7","A19","A31","A43"],
    "gal2": ["A8","A20","A32","A44"],
    "gal3": ["A9","A21","A33","A45"],
    "gal4": ["A10","A22","A34","A46"],
    "gal5": ["A11","A23","A35","A47"],
    "gal6": ["A12","A24","A36","A48"],
}

results = []

for label, cols in sets.items():
    # Merge all genes across the 4 timepoints
    merged = pd.concat([df[c].dropna() for c in cols])

    total_genes = len(merged)                     # total entries (with duplicates)
    unique_genes = merged.nunique()               # number of unique genes
    percent_genes = (unique_genes / 6281) * 100   # percentage

    row = {
        "Set": label,
        "Total_genes": total_genes,
        "Unique_genes": unique_genes,
        "Percent_genes": round(percent_genes, 2)
    }
    results.append(row)

# Save output
out_df = pd.DataFrame(results)
out_df.to_csv(output_file, sep="\t", index=False)

print(f"Results saved to {output_file}")


