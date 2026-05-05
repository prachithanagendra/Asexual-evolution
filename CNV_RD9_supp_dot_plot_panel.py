#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os

# === Input file paths ===
file_300 = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/3_merged_calls/merged_CNV_state_glu_300.tsv"
file_600 = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/3_merged_calls/merged_CNV_state_glu_600.tsv"
file_900 = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/3_merged_calls/merged_CNV_state_glu_900.tsv"
file_1200 = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/3_merged_calls/merged_CNV_state_glu_1200.tsv"

# === Chromosome order (yeast 16 chromosomes) ===
chrom_order = [
    'I','II','III','IV','V','VI','VII','VIII',
    'IX','X','XI','XII','XIII','XIV','XV','XVI'
]

# Color map
chrom_to_color = {chrom: cm.tab20(i % 20) for i, chrom in enumerate(chrom_order)}


# ==========================================================
# FILE PROCESSING
# ==========================================================
def process_file(filepath, col_name):
    df = pd.read_csv(filepath, sep="\t")

    if col_name not in df.columns:
        raise ValueError(
            f"Column '{col_name}' not found in {filepath}. "
            f"Available columns: {list(df.columns)}"
        )

    df['value'] = df[col_name]

    # Standardize chromosome names
    df['chrom_id'] = df['chrom'].astype(str).str.replace('chr', '', regex=False)

    return df[['chrom_id', 'window', 'value']]


# ==========================================================
# GLOBAL X POSITION MAPPING
# ==========================================================
def assign_global_xpos(df_base):
    xpos_map = {}
    boundaries = {}
    current_offset = 0

    for i, chrom in enumerate(chrom_order):
        chrom_df = df_base[df_base['chrom_id'] == chrom].sort_values('window')

        if chrom_df.empty:
            print(f"Warning: No data found for chromosome {chrom}")
            continue

        chrom_len = chrom_df['window'].max()
        boundaries[chrom] = (current_offset, current_offset + chrom_len)

        for _, row in chrom_df.iterrows():
            xpos_map[(chrom, row['window'])] = current_offset + row['window']

        if i < len(chrom_order) - 1:
            current_offset += chrom_len + 50000
        else:
            current_offset += chrom_len

    return xpos_map, boundaries


def apply_xpos(df, xpos_map):
    df = df.copy()

    df['x_pos'] = df.apply(
        lambda row: xpos_map.get((row['chrom_id'], row['window']), np.nan),
        axis=1
    )

    # Assign colors safely
    df['color'] = df['chrom_id'].map(chrom_to_color)

    # Remove rows with missing positions
    df = df.dropna(subset=['x_pos'])

    # Replace unmapped chromosomes with black (prevents crash)
    df['color'] = df['color'].fillna('black')

    return df


# ==========================================================
# INDIVIDUAL REPLICATE PLOTS
# ==========================================================
def plot_replicate(rep_name, col_map, output_dir, xpos_map, chrom_boundaries):

    df_300 = process_file(file_300, col_map["300"])
    df_600 = process_file(file_600, col_map["600"])
    df_900 = process_file(file_900, col_map["900"])
    df_1200 = process_file(file_1200, col_map["1200"])

    dfs = {
        "300 generations": apply_xpos(df_300, xpos_map),
        "600 generations": apply_xpos(df_600, xpos_map),
        "900 generations": apply_xpos(df_900, xpos_map),
        "1200 generations": apply_xpos(df_1200, xpos_map),
    }

    ymin, ymax = 0.4, 1.8
    fig, axes = plt.subplots(4, 1, figsize=(18, 12), sharex=True)

    for ax, (label, df) in zip(axes, dfs.items()):
        ax.scatter(df['x_pos'], df['value'], color=df['color'], s=8)
        ax.set_title(f"{rep_name} - {label}", fontsize=14)
        ax.set_ylim(ymin, ymax)

        for chrom in chrom_order[:-1]:
            if chrom in chrom_boundaries:
                _, end = chrom_boundaries[chrom]
                ax.axvline(x=end + 25000, color='black', linewidth=0.8)

    # X-axis ticks
    xticks, xticklabels = [], []
    for chrom in chrom_order:
        if chrom in chrom_boundaries:
            start, end = chrom_boundaries[chrom]
            xticks.append((start + end) // 2)
            xticklabels.append(chrom)

    axes[-1].set_xticks(xticks)
    axes[-1].set_xticklabels(xticklabels, fontsize=12)
    axes[-1].set_xlabel("Chromosome")

    plt.tight_layout()
    outfile = os.path.join(output_dir, f"{rep_name}_CNV_state_plot.png")
    plt.savefig(outfile, dpi=300)
    plt.close(fig)

    return dfs


# ==========================================================
# MAIN
# ==========================================================
output_dir = "/home1/NGF_SK1_reference/NGF_SK1_ref_analysis/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/4_avg_CNV_plots/glu"
os.makedirs(output_dir, exist_ok=True)

replicate_columns = {
    "glu1": {"300": "A1", "600": "A13", "900": "A25", "1200": "A37"},
    "glu2": {"300": "A2", "600": "A14", "900": "A26", "1200": "A38"},
    "glu3": {"300": "A3", "600": "A15", "900": "A27", "1200": "A39"},
    "glu4": {"300": "A4", "600": "A16", "900": "A28", "1200": "A40"},
    "glu5": {"300": "A5", "600": "A17", "900": "A29", "1200": "A41"},
    "glu6": {"300": "A6", "600": "A18", "900": "A30", "1200": "A42"},
}

# Use first replicate to build chromosome layout
df_base = process_file(file_300, "A1")
xpos_map, chrom_boundaries = assign_global_xpos(df_base)

all_data = {}
for rep_name, col_map in replicate_columns.items():
    all_data[rep_name] = plot_replicate(
        rep_name, col_map, output_dir, xpos_map, chrom_boundaries
    )

print(f"Saved individual replicate plots in {output_dir}")


# ==========================================================
# BIG PANEL (6 x 4)
# ==========================================================
time_labels = ["300 generations", "600 generations",
               "900 generations", "1200 generations"]

ymin, ymax = 0.4, 1.8
fig, axes = plt.subplots(6, 4, figsize=(36, 24),
                         sharex=True, sharey=True)

for i, (rep_name, dfs) in enumerate(all_data.items()):
    for j, label in enumerate(time_labels):
        ax = axes[i, j]
        df = dfs[label]
        ax.scatter(df['x_pos'], df['value'], color=df['color'], s=8)

        if i == 0:
            ax.set_title(label, fontsize=18)
        if j == 0:
            ax.set_ylabel(rep_name, fontsize=18)

        ax.set_ylim(ymin, ymax)

        for chrom in chrom_order[:-1]:
            if chrom in chrom_boundaries:
                _, end = chrom_boundaries[chrom]
                ax.axvline(x=end + 25000, color='black', linewidth=0.8)

# Bottom xticks only
xticks, xticklabels = [], []
for chrom in chrom_order:
    if chrom in chrom_boundaries:
        start, end = chrom_boundaries[chrom]
        xticks.append((start + end) // 2)
        xticklabels.append(chrom)

for ax in axes[-1, :]:
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=20)

fig.text(0.5, 0.04, "Chromosome", ha='center', fontsize=20)
fig.text(0.04, 0.5, "CNV state", va='center',
         rotation='vertical', fontsize=20)

plt.tight_layout(rect=[0.06, 0.06, 1, 0.97])

outfile_png = os.path.join(output_dir,
                           "glu_allreplicates_CNV_state_bigpanel.png")
outfile_pdf = os.path.join(output_dir,
                           "glu_allreplicates_CNV_state_bigpanel.pdf")

plt.savefig(outfile_png, dpi=300)
plt.savefig(outfile_pdf)
plt.close(fig)

print(f"Saved big panel plot at {outfile_png} and {outfile_pdf}")

