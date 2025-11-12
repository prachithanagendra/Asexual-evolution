import os
import subprocess
import argparse

def run_command(command, step_folder):
    """Runs a shell command and stores output in the respective folder."""
    print(f"Running: {command}")
    subprocess.run(command, shell=True, check=True, cwd=step_folder)

def main():
    parser = argparse.ArgumentParser(description="WGS Pipeline")
    parser.add_argument("--base_dir", required=True, help="Path to input sequencing data")
    parser.add_argument("--output_base", required=True, help="Path to output directory")
    parser.add_argument("--reference_genome", required=True, help="Path to reference genome")
    args = parser.parse_args()

    base_dir = args.base_dir
    output_base = args.output_base
    reference_genome = args.reference_genome
    picard_cmd = "java -jar /home/prachitha/WGS/picard.jar"
    gatk_cmd = "java -jar /home/prachitha/WGS/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar"

    samples = [f"A{i}" for i in range(1, 49)]

    folder_names = {
        1: "1_Trimming",
        2: "2_Quality_Check",
        3: "3_Adapter_Removal",
        4: "4_Alignment",
        5: "5_Filtering",
        6: "6_Sorting",
        7: "7_Read_Group_And_MarkDuplicates",
        8: "8_Sorting_And_Indexing",
        9: "9_Variant_Calling",
        10: "10_Assembled_Reads"
    }

    for step in folder_names:
        step_folder = os.path.join(output_base, folder_names[step])
        os.makedirs(step_folder, exist_ok=True)

    for sample in samples:
        r1_file = os.path.join(base_dir, f"{sample}_WGS_BatchEXT24_R1.fastq.gz")
        r2_file = os.path.join(base_dir, f"{sample}_WGS_BatchEXT24_R2.fastq.gz")

        # Step 1: Fastp trimming
        trimmed_r1 = os.path.join(output_base, folder_names[1], f"{sample}_R1_trimmed.fastq")
        trimmed_r2 = os.path.join(output_base, folder_names[1], f"{sample}_R2_trimmed.fastq")
        command = f"fastp -i {r1_file} -I {r2_file} -o {trimmed_r1} -O {trimmed_r2} --html {sample}_fastp_report.html --json {sample}_fastp_report.json --thread 4"
        run_command(command, os.path.join(output_base, folder_names[1]))

        # Step 2: FastQC
        command = f"fastqc {trimmed_r1} {trimmed_r2}"
        run_command(command, os.path.join(output_base, folder_names[2]))

        # Step 3: Adapter removal
        nextera_trimmed_r1 = os.path.join(output_base, folder_names[3], f"{sample}_nextera_trimmed_R1.fastq.gz")
        nextera_trimmed_r2 = os.path.join(output_base, folder_names[3], f"{sample}_nextera_trimmed_R2.fastq.gz")
        command = f"fastp -i {trimmed_r1} -I {trimmed_r2} -o {nextera_trimmed_r1} -O {nextera_trimmed_r2} --detect_adapter_for_pe -q 20 -l 36"
        run_command(command, os.path.join(output_base, folder_names[3]))

        # Step 4: FastQC after adapter removal
        command = f"fastqc {nextera_trimmed_r1} {nextera_trimmed_r2}"
        run_command(command, os.path.join(output_base, folder_names[2]))

        # Step 5: Alignment using BWA
        raw_sam = os.path.join(output_base, folder_names[4], f"{sample}_raw.sam")
        command = f"bwa mem -t 4 -k 32 -M {reference_genome} {nextera_trimmed_r1} {nextera_trimmed_r2} > {raw_sam}"
        run_command(command, os.path.join(output_base, folder_names[4]))

        # Step 6: Convert SAM to BAM
        raw_bam = os.path.join(output_base, folder_names[5], f"{sample}_raw.bam")
        command = f"samtools view -bS {raw_sam} > {raw_bam}"
        run_command(command, os.path.join(output_base, folder_names[5]))

        # Step 7: Sorting BAM file
        sorted_bam = os.path.join(output_base, folder_names[6], f"{sample}_sorted.bam")
        command = f"samtools sort {raw_bam} -o {sorted_bam}"
        run_command(command, os.path.join(output_base, folder_names[6]))

        # Step 8: Add Read Groups and Mark Duplicates
        rg_bam = os.path.join(output_base, folder_names[7], f"{sample}.rg.bam")
        command = f"{picard_cmd} AddOrReplaceReadGroups -I {sorted_bam} -O {rg_bam} -ID {sample} -LB {sample} -PL ILLUMINA -SM {sample} -PU {sample}"
        run_command(command, os.path.join(output_base, folder_names[7]))

        marked_dup_bam = os.path.join(output_base, folder_names[7], f"{sample}.marked_dup.bam")
        metrics_file = os.path.join(output_base, folder_names[7], f"{sample}_metrics_duplicate.txt")
        command = f"{picard_cmd} MarkDuplicates -I {rg_bam} -O {marked_dup_bam} -M {metrics_file}"
        run_command(command, os.path.join(output_base, folder_names[7]))

        # Step 9: Sorting and Indexing BAM file
        sorted1_bam = os.path.join(output_base, folder_names[8], f"{sample}_sorted1.bam")
        command = f"{picard_cmd} SortSam -I {marked_dup_bam} -O {sorted1_bam} -SO coordinate"
        run_command(command, os.path.join(output_base, folder_names[8]))

        command = f"samtools index {sorted1_bam}"
        run_command(command, os.path.join(output_base, folder_names[8]))

        # Step 10: Variant Calling using GATK
        vcf_raw = os.path.join(output_base, folder_names[9], f"{sample}_GATK.raw.vcf")
        assembled_reads_bam = os.path.join(output_base, folder_names[10], f"{sample}_assembled_reads.bam")
        assembled_regions_txt = os.path.join(output_base, folder_names[10], f"{sample}_assembled_regions.txt")
        command = f"{gatk_cmd} HaplotypeCaller --output {vcf_raw} --input {sorted1_bam} --reference {reference_genome} --ploidy 2 --bam-output {assembled_reads_bam} --assembly-region-out {assembled_regions_txt}"
        run_command(command, os.path.join(output_base, folder_names[9]))

if __name__ == "__main__":
    main()

