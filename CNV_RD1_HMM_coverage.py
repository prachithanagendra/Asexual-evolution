import os
import subprocess
import pandas as pd
import numpy as np
from hmmlearn.hmm import GaussianHMM

# === SETTINGS ===
bam_dir = "/8_Sorting_And_Indexing_A"
output_dir = "/25_Final_CNV_analysis/A_Read_depth"
window_size = 500
ancestor_id = "A0_sorted1"

os.makedirs(f"{output_dir}/1_depth", exist_ok=True)
os.makedirs(f"{output_dir}/2_cnv_calls", exist_ok=True)

# === STEP 1: CALCULATE DEPTH WITH SAMTOOLS ===
def generate_depth_file(bam_path, output_path):
    cmd = f"samtools depth {bam_path} > {output_path}"
    subprocess.run(cmd, shell=True, check=True)

# === STEP 2: BIN DEPTH INTO 500bp WINDOWS ===
def bin_depth_file(depth_file, window_size):
    df = pd.read_csv(depth_file, sep="\t", names=["chrom", "pos", "depth"])
    df["window"] = (df["pos"] // window_size) * window_size
    binned = df.groupby(["chrom", "window"])["depth"].median().reset_index()
    binned["end"] = binned["window"] + window_size
    return binned[["chrom", "window", "end", "depth"]]

# === STEP 3: NORMALIZE DEPTH ===
def compute_relative_depth(df):
    median = df["depth"].median()
    df["relative_depth"] = df["depth"] / median
    return df

# === STEP 4: STANDARDIZE USING ANCESTOR RELATIVE DEPTH ===
def standardize_by_ancestor(sample_df, ancestor_df):
    merged = sample_df.merge(ancestor_df, on=["chrom", "window", "end"], suffixes=("", "_ancestor"))
    merged = merged[merged["relative_depth_ancestor"] >= 0.25].copy()
    merged["standardized"] = merged["relative_depth"] / merged["relative_depth_ancestor"]
    return merged

# === STEP 5: HMM CNV DETECTION ===
def run_hmm(df, sample_name, out_path):
    states = [0, 0.5, 1, 1.5, 2, 3, 4]
    X = df["standardized"].values.reshape(-1, 1)
    baseline_var = np.var(X)
    means = np.array([[s] for s in states])
    covars = np.array([[baseline_var * 0.5 if s == 0 else baseline_var * s] for s in states])

    model = GaussianHMM(n_components=len(states), covariance_type="diag", init_params="")
    model.means_ = means
    model.covars_ = covars
    model.startprob_ = np.array([0.01] * len(states))
    model.startprob_[states.index(1)] = 0.94
    transmat = np.full((len(states), len(states)), 0.0001)
    np.fill_diagonal(transmat, 0.9994)
    model.transmat_ = transmat

    df["CNV_state"] = [states[s] for s in model.predict(X)]
    df.to_csv(out_path, sep="\t", index=False)

# === MAIN WORKFLOW ===
binned_depths = {}

# Process ancestor BAM first
bam_files = [f for f in os.listdir(bam_dir) if f.endswith(".bam")]
ancestor_bam = f"{ancestor_id}.bam"
if ancestor_bam not in bam_files:
    raise ValueError(f"Ancestor BAM file {ancestor_bam} not found in {bam_dir}")

bam_files.remove(ancestor_bam)
bam_files.insert(0, ancestor_bam)

for bam in bam_files:
    sample = bam.replace(".bam", "")
    bam_path = os.path.join(bam_dir, bam)
    depth_path = f"{output_dir}/1_depth/{sample}.depth"

    if os.path.exists(depth_path):
        print(f"Depth file already exists for {sample}, skipping samtools depth.")
    else:
        print(f"Running samtools depth for {sample}...")
        generate_depth_file(bam_path, depth_path)

    binned = bin_depth_file(depth_path, window_size)
    binned = compute_relative_depth(binned)
    binned_depths[sample] = binned

# Get ancestor baseline
ancestor_df = binned_depths[ancestor_id].copy()
ancestor_df = ancestor_df[["chrom", "window", "end", "relative_depth"]]
ancestor_df = ancestor_df.rename(columns={"relative_depth": "relative_depth_ancestor"})

# Process all samples except ancestor
for sample, df in binned_depths.items():
    if sample == ancestor_id:
        continue
    standardized_df = standardize_by_ancestor(df, ancestor_df)
    out_path = f"{output_dir}/2_cnv_calls/{sample}_cnv_calls.tsv"
    run_hmm(standardized_df, sample, out_path)
    print(f"✔ CNV calls saved for {sample} → {out_path}")

