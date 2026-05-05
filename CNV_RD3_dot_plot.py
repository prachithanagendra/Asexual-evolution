import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# === Input file paths ===
file_300 = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/3_merged_calls/merged_CNV_state_gal_300.tsv"
file_600 = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/3_merged_calls/merged_CNV_state_gal_600.tsv"
file_900 = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/3_merged_calls/merged_CNV_state_gal_900.tsv"
file_1200 = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/3_merged_calls/merged_CNV_state_gal_1200.tsv"

# === Chromosome label mapping (Roman numerals) ===
chrom_label_map = {
    'I': 'I', 'II': 'II', 'III': 'III', 'IV': 'IV',
    'V': 'V', 'VI': 'VI', 'VII': 'VII', 'VIII': 'VIII',
    'IX': 'IX', 'X': 'X', 'XI': 'XI', 'XII': 'XII',
    'XIII': 'XIII', 'XIV': 'XIV', 'XV': 'XV', 'XVI': 'XVI'
}

# === Preprocess file ===
def process_file(filepath):
    df = pd.read_csv(filepath, sep="\t")

    sample_cols = df.columns[-6:]
    df['avg'] = df[sample_cols].mean(axis=1)

    # FIX: chromosomes are chrI, chrII → strip "chr"
    df['chrom_id'] = df['chrom'].astype(str).str.replace("^chr", "", regex=True)
    df['chrom_label'] = df['chrom_id']

    return df[['chrom_id', 'chrom_label', 'window', 'avg']]

# === Load all datasets ===
df_300 = process_file(file_300)
df_600 = process_file(file_600)
df_900 = process_file(file_900)
df_1200 = process_file(file_1200)

# === Define chromosome order and color mapping ===
chrom_order = list(chrom_label_map.keys())
chrom_to_color = {chrom: cm.tab20(i % 20) for i, chrom in enumerate(chrom_order)}


# === Assign consistent x_pos and chromosome boundaries ===
def assign_global_xpos(df_base):
    xpos_map = {}
    boundaries = {}
    current_offset = 0

    for i, chrom in enumerate(chrom_order):
        chrom_df = df_base[df_base['chrom_id'] == chrom].sort_values('window')
        chrom_len = chrom_df['window'].max()

        boundaries[chrom] = (current_offset, current_offset + chrom_len)

        for _, row in chrom_df.iterrows():
            xpos_map[(chrom, row['window'])] = current_offset + row['window']

        if i < len(chrom_order) - 1:
            current_offset += chrom_len + 50000
        else:
            current_offset += chrom_len

    return xpos_map, boundaries

xpos_map, chrom_boundaries = assign_global_xpos(df_300)

# === Apply consistent x_pos to all datasets ===
def apply_xpos(df, xpos_map):
    df = df.copy()
    df['x_pos'] = df.apply(lambda row: xpos_map.get((row['chrom_id'], row['window']), np.nan), axis=1)
    df['color'] = df['chrom_id'].map(chrom_to_color)
    return df

df_300 = apply_xpos(df_300, xpos_map)
df_600 = apply_xpos(df_600, xpos_map)
df_900 = apply_xpos(df_900, xpos_map)
df_1200 = apply_xpos(df_1200, xpos_map)

# === Fixed y-axis limits ===
ymin, ymax = 0.4, 1.8

# === Create subplots ===
fig, axes = plt.subplots(4, 1, figsize=(18, 12), sharex=True)
fig.text(0.04, 0.5, 'CNV state', va='center', rotation='vertical', fontsize=14, fontname='Arial')

# === Timepoint mapping ===
timepoints = {
    '300 generations': df_300,
    '600 generations': df_600,
    '900 generations': df_900,
    '1200 generations': df_1200
}

# === Plot each panel ===
for ax, (label, df) in zip(axes, timepoints.items()):
    ax.scatter(df['x_pos'], df['avg'], color=df['color'], s=8)
    ax.set_title(label, fontsize=14)
    ax.set_ylim(ymin, ymax)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_yticks(np.arange(ymin, ymax + 0.01, 0.2))

    for chrom in chrom_order[:-1]:
        _, end = chrom_boundaries[chrom]
        ax.axvline(x=end + 25000, color='black', linestyle='-', linewidth=0.8)

# === Set chromosome labels on bottom axis ===
xticks = []
xticklabels = []
for chrom in chrom_order:
    start, end = chrom_boundaries[chrom]
    xticks.append((start + end) // 2)
    xticklabels.append(chrom_label_map[chrom])

axes[-1].set_xticks(xticks)
axes[-1].set_xticklabels(xticklabels, fontsize=12, fontname='Arial')
axes[-1].set_xlabel("Chromosome", fontsize=14, fontname='Arial')

start_x = chrom_boundaries[chrom_order[0]][0]
end_x = chrom_boundaries[chrom_order[-1]][1]
for ax in axes:
    ax.set_xlim(start_x, end_x)
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')


# === Save and show ===
plt.tight_layout(rect=[0.06, 0.03, 1, 0.97])
plt.savefig(
    "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/4_avg_CNV_plots/galactose_CNV_state_plot.png",
    dpi=300,
    bbox_inches='tight'
)
plt.show()

