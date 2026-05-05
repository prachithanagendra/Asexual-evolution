#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
from pathlib import Path

# === USER CONFIG ===
input_folder = Path("/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/5_CNV_separated/duplication")
gff_file = Path("/SK1_lifted_annotation.gff3")
output_base_dir = Path("/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/6_CNV_separated_annotated/duplication_annotated")

# Create output directories
output_base_dir.mkdir(parents=True, exist_ok=True)
unique_genes_dir = output_base_dir / "unique_genes_dup"
unique_genes_dir.mkdir(parents=True, exist_ok=True)

# ==========================================================
# STEP 0: Extract genes.bed (only once)
# ==========================================================

genes_bed = output_base_dir / "genes.bed"

if not genes_bed.exists():
    print(" Extracting genes from GFF file...")
    with open(gff_file) as fin, open(genes_bed, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand, _, attr = parts
            attrs = {k: v for k, v in (item.split('=') for item in attr.split(';') if '=' in item)}
            gene_id = attrs.get('ID', '.')
            gene_name = attrs.get('Name', '.')
            fout.write(f"{chrom}\t{int(start)-1}\t{end}\t{gene_id}\t{gene_name}\t{strand}\n")

# ==========================================================
# LOOP OVER A1–A48
# ==========================================================

for i in range(1, 49):

    sample_name = f"A{i}"
    input_filename = f"{sample_name}_sorted1_cnv_calls.tsv"
    window_tsv = input_folder / input_filename

    if not window_tsv.exists():
        print(f" Missing file, skipping: {window_tsv}")
        continue

    print(f"\n Processing {sample_name} ...")

    # Define output files
    windows_bed = output_base_dir / f"{sample_name}_windows.bed"
    merged_bed = output_base_dir / f"{sample_name}_merged_blocks.bed"
    intersect_file = output_base_dir / f"{sample_name}_full_gene_overlap.tsv"
    parsed_output = output_base_dir / f"{sample_name}_parsed_annotation.tsv"
    filtered_output = output_base_dir / f"{sample_name}_parsed_annotation_filtered.tsv"
    unique_genes_txt = unique_genes_dir / f"{sample_name}_unique_gene_names_dup.txt"

    # ==========================================================
    # STEP 1: Convert TSV → BED
    # ==========================================================

    with open(window_tsv) as fin, open(windows_bed, "w") as fout:
        next(fin)  # skip header
        for line in fin:
            chrom, start, end = line.strip().split('\t')[:3]
            if chrom == "BK006947.2":
                chrom = "BK006947.3"
            start_bed = max(0, int(start) - 1)
            fout.write(f"{chrom}\t{start_bed}\t{end}\n")

    # ==========================================================
    # STEP 2: Merge adjacent windows
    # ==========================================================

    subprocess.run([
        "bedtools", "merge",
        "-i", str(windows_bed)
    ], stdout=open(merged_bed, "w"), check=True)

    # ==========================================================
    # STEP 3: Intersect genes fully contained in CNV blocks
    # ==========================================================

    subprocess.run([
        "bedtools", "intersect",
        "-a", str(genes_bed),
        "-b", str(merged_bed),
        "-wa", "-u",
        "-f", "1.0"
    ], stdout=open(intersect_file, "w"), check=True)

    # ==========================================================
    # STEP 4: Parse intersect results
    # ==========================================================

    if os.path.getsize(intersect_file) == 0:
        print(f"⚠ No fully duplicated genes found for {sample_name}")
        continue

    cols = [
        "gene_chrom", "gene_start", "gene_end",
        "gene_id", "gene_name", "strand"
    ]

    df = pd.read_csv(intersect_file, sep='\t', header=None, names=cols)
    df.drop_duplicates(inplace=True)
    df.to_csv(parsed_output, sep='\t', index=False)

    # ==========================================================
    # STEP 5: Filter gene_id != gene_name
    # ==========================================================

    df["core_gene_id"] = df["gene_id"].str.replace(
        r"^(gene-|transcript-|.*?:)?", "", regex=True
    )

    filtered_df = df[df["core_gene_id"] != df["gene_name"]].copy()
    filtered_df.drop(columns=["core_gene_id"], inplace=True)
    filtered_df.to_csv(filtered_output, sep='\t', index=False)

    # ==========================================================
    # STEP 6: Extract unique gene names
    # ==========================================================

    unique_genes = sorted(filtered_df["gene_name"].unique())

    with open(unique_genes_txt, "w") as fout:
        fout.write("\n".join(unique_genes))

    print(f" Completed {sample_name}")
    print(f"  Fully duplicated genes: {len(unique_genes)}")

print("\n All samples processed successfully.")

