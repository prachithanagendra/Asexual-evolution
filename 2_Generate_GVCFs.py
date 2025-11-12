import os
import subprocess

# User-defined GATK command
gatk_cmd = "java -jar /home/prachitha/WGS/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"

# Prompt user for inputs
bam_dir = input("Enter the path to the folder containing BAM files: ").strip()
output_dir = input("Enter the path to the output folder for GVCFs: ").strip()
reference = input("Enter the path to the reference genome (FASTA): ").strip()

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Loop over sample numbers
for i in range(1, 49):
    sample_id = f"A{i}"
    bam_file = os.path.join(bam_dir, f"{sample_id}_sorted1.bam")
    output_file = os.path.join(output_dir, f"{sample_id}.g.vcf")

    if not os.path.isfile(bam_file):
        print(f"⚠ Skipping {sample_id}: BAM file not found at {bam_file}")
        continue

    print(f"Processing {sample_id}...")

    cmd = f'{gatk_cmd} HaplotypeCaller -R "{reference}" -I "{bam_file}" -O "{output_file}" -ERC GVCF'

    try:
        subprocess.run(cmd, check=True, shell=True)
        print(f"✔ Finished: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error processing {sample_id}: {e}")

print("✅ All available BAM files processed.")

