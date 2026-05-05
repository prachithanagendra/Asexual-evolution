import pandas as pd
from pathlib import Path

# User configurable paths
input_folder = Path("/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/2_cnv_calls")
output_folder = Path("/25_Final_CNV_analysis/A_Read_depth/8_redone_10032026/5_CNV_separated")
deletion_folder = output_folder / "deletion"
duplication_folder = output_folder / "duplication"

# Create output subfolders if not exist
deletion_folder.mkdir(parents=True, exist_ok=True)
duplication_folder.mkdir(parents=True, exist_ok=True)

for i in range(1, 49):
    filename = f"A{i}_sorted1_cnv_calls.tsv"
    input_file = input_folder / filename
    
    if not input_file.exists():
        print(f"Warning: {filename} not found, skipping.")
        continue
    
    # Read TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    # Filter rows based on CNV_state
    deletion_df = df[df['CNV_state'] == 0.5]
    duplication_df = df[df['CNV_state'] == 1.5]
    
    # Write to respective output files
    deletion_output_file = deletion_folder / filename
    duplication_output_file = duplication_folder / filename
    
    deletion_df.to_csv(deletion_output_file, sep='\t', index=False)
    duplication_df.to_csv(duplication_output_file, sep='\t', index=False)
    
    print(f"Processed {filename}: {len(deletion_df)} deletions, {len(duplication_df)} duplications")

print("All files processed.")

