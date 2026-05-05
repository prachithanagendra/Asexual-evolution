import os
import pandas as pd
import numpy as np
from hmmlearn.hmm import GaussianHMM

# =========================================================
# SETTINGS
# =========================================================
depth_dir = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/1_depth"
output_dir = "/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/except_telomere"
window_size = 500
ancestor_id = "A0_sorted1"

os.makedirs(f"{output_dir}/2_cnv_calls", exist_ok=True)

# =========================================================
# STEP 1: EXCLUDE FIRST/LAST 500 bp AND BIN INTO WINDOWS
# =========================================================
def bin_depth_file(depth_file, window_size):
    df = pd.read_csv(
        depth_file,
        sep="\t",
        names=["chrom", "pos", "depth"]
    )

    filtered_chromosomes = []

    for chrom, chrom_df in df.groupby("chrom"):
        max_pos = chrom_df["pos"].max()

        # Exclude first 500 bp and last 500 bp
        chrom_df = chrom_df[
            (chrom_df["pos"] > window_size) &
            (chrom_df["pos"] <= (max_pos - window_size))
        ].copy()

        filtered_chromosomes.append(chrom_df)

    df = pd.concat(filtered_chromosomes, ignore_index=True)

    # Create 500 bp windows
    df["window"] = (df["pos"] // window_size) * window_size

    # Median depth per window
    binned = (
        df.groupby(["chrom", "window"])["depth"]
        .median()
        .reset_index()
    )

    binned["end"] = binned["window"] + window_size

    return binned[["chrom", "window", "end", "depth"]]

# =========================================================
# STEP 2: COMPUTE RELATIVE DEPTH
# =========================================================
def compute_relative_depth(df):
    median_depth = df["depth"].median()
    df["relative_depth"] = df["depth"] / median_depth
    return df

# =========================================================
# STEP 3: STANDARDIZE USING ANCESTOR
# =========================================================
def standardize_by_ancestor(sample_df, ancestor_df):
    merged = sample_df.merge(
        ancestor_df,
        on=["chrom", "window", "end"],
        suffixes=("", "_ancestor")
    )

    merged = merged[
        merged["relative_depth_ancestor"] >= 0.25
    ].copy()

    merged["standardized"] = (
        merged["relative_depth"] /
        merged["relative_depth_ancestor"]
    )

    return merged

# =========================================================
# STEP 4: HMM CNV DETECTION
# =========================================================
def run_hmm(df, out_path):
    states = [0, 0.5, 1, 1.5, 2, 3, 4]

    X = df["standardized"].values.reshape(-1, 1)

    baseline_var = np.var(X)

    means = np.array([[s] for s in states])

    covars = np.array([
        [baseline_var * 0.5 if s == 0 else baseline_var * s]
        for s in states
    ])

    model = GaussianHMM(
        n_components=len(states),
        covariance_type="diag",
        init_params=""
    )

    model.means_ = means
    model.covars_ = covars

    # Start probabilities
    model.startprob_ = np.array([0.01] * len(states))
    model.startprob_[states.index(1)] = 0.94

    # Transition matrix
    transmat = np.full(
        (len(states), len(states)),
        0.0001
    )
    np.fill_diagonal(transmat, 0.9994)

    model.transmat_ = transmat

    predicted_states = model.predict(X)

    df["CNV_state"] = [
        states[s] for s in predicted_states
    ]

    df.to_csv(out_path, sep="\t", index=False)

# =========================================================
# MAIN WORKFLOW
# =========================================================
binned_depths = {}

# Get all depth files
depth_files = [
    f for f in os.listdir(depth_dir)
    if f.endswith(".depth")
]

# Ensure ancestor is processed first
ancestor_depth = f"{ancestor_id}.depth"

if ancestor_depth not in depth_files:
    raise ValueError(
        f"Ancestor depth file {ancestor_depth} not found in {depth_dir}"
    )

depth_files.remove(ancestor_depth)
depth_files.insert(0, ancestor_depth)

# =========================================================
# PROCESS ALL DEPTH FILES
# =========================================================
for depth_file in depth_files:
    sample = depth_file.replace(".depth", "")
    depth_path = os.path.join(depth_dir, depth_file)

    print(f"Processing {sample}...")

    binned = bin_depth_file(depth_path, window_size)
    binned = compute_relative_depth(binned)

    binned_depths[sample] = binned

# =========================================================
# PREPARE ANCESTOR BASELINE
# =========================================================
ancestor_df = binned_depths[ancestor_id].copy()

ancestor_df = ancestor_df[
    ["chrom", "window", "end", "relative_depth"]
]

ancestor_df = ancestor_df.rename(
    columns={
        "relative_depth": "relative_depth_ancestor"
    }
)

# =========================================================
# RUN CNV CALLING
# =========================================================
for sample, df in binned_depths.items():
    if sample == ancestor_id:
        continue

    print(f"Running HMM for {sample}...")

    standardized_df = standardize_by_ancestor(
        df,
        ancestor_df
    )

    out_path = (
        f"{output_dir}/2_cnv_calls/"
        f"{sample}_cnv_calls.tsv"
    )

    run_hmm(standardized_df, out_path)

    print(f"✔ CNV calls saved → {out_path}")
