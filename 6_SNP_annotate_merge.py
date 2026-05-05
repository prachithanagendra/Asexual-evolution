#!/usr/bin/env python3

import os
import pandas as pd
import subprocess

# ==========================================================
# USER INPUTS
# ==========================================================

INPUT_CSV      = "/14_SNP_analysis/3_fixed_variants/tsv/gal6_B.csv"
JOINT_VCF      = "/2_Joint_variantcalling/joint_genotyped_all_samples.vcf"
AF_TSV         = "/14_SNP_analysis/3_fixed_variants/tsv/B12_B24.tsv"

WORK_FOLDER    = "/14_SNP_analysis/4_fixed_annotated"
FINAL_OUTPUT   = "/14_SNP_analysis/5_fixed_annotated_merged"

os.makedirs(WORK_FOLDER, exist_ok=True)
os.makedirs(FINAL_OUTPUT, exist_ok=True)

# intermediate files
FILTERED_VCF   = os.path.join(WORK_FOLDER, "gal6_B_filtered.vcf")
ANNOTATED_VCF  = os.path.join(WORK_FOLDER, "gal6_B_annotated.vcf")
ANNOTATED_TSV  = os.path.join(WORK_FOLDER, "gal6_B_annotated.tsv")
FILTERED_ANN   = os.path.join(WORK_FOLDER, "gal6_B_filtered_annotation.tsv")

FINAL_TSV      = os.path.join(FINAL_OUTPUT, "gal6_B_annotation_AF_merged.tsv")


# ==========================================================
# STEP 1: CSV → extract matching variants from joint VCF
# ==========================================================

def csv_to_vcf(csv_file, joint_vcf, output_vcf):

    df = pd.read_csv(csv_file, sep="\t")

    variant_set = set(zip(df["CHROM"], df["POS"], df["REF"], df["ALT"]))

    filtered_lines = []

    with open(joint_vcf) as vcf:
        for line in vcf:

            if line.startswith("#"):
                filtered_lines.append(line)
                continue

            chrom, pos, _, ref, alt = line.strip().split("\t")[:5]

            if (chrom, int(pos), ref, alt) in variant_set:
                filtered_lines.append(line)

    with open(output_vcf, "w") as out:
        out.writelines(filtered_lines)


# ==========================================================
# STEP 2: run snpEff
# ==========================================================

def run_snpeff(input_vcf, output_vcf):

    print("Running snpEff annotation")

    cmd = f"java -jar snpEff.jar ann SK1_custom {input_vcf} > {output_vcf}"

    subprocess.run(cmd, shell=True, check=True)


# ==========================================================
# STEP 3: annotated VCF → TSV
# ==========================================================

def vcf_to_tsv(input_vcf, output_tsv):

    awk_cmd = r"""awk 'BEGIN {OFS="\t"; print "CHROM","POS","REF","ALT","QUAL","Allele","Consequence","Impact","GeneName","GeneID","FeatureType","FeatureID","TranscriptBioType"} {n = split($6, annots, ","); for(i=1;i<=n;i++){split(annots[i], a, "|"); print $1,$2,$3,$4,$5,a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8]}}'"""

    cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/ANN\n' {input_vcf} | {awk_cmd} > {output_tsv}"

    subprocess.run(cmd, shell=True, check=True)


# ==========================================================
# STEP 4: keep all annotations ONLY if first annotation is
# upstream_gene_variant or downstream_gene_variant
# ==========================================================

def keep_first_annotation(input_tsv, output_tsv):

    df = pd.read_csv(input_tsv, sep="\t")

    def filter_annotations(group):

        first_cons = str(group.iloc[0]["Consequence"])

        if ("upstream_gene_variant" in first_cons) or ("downstream_gene_variant" in first_cons):
            return group

        # otherwise keep only first annotation
        return group.iloc[[0]]

    filtered_df = (
        df.groupby(["CHROM","POS","REF","ALT"], group_keys=False)
        .apply(filter_annotations)
    )

    filtered_df.to_csv(output_tsv, sep="\t", index=False)

# ==========================================================
# STEP 5: merge with allele frequency table
# ==========================================================

def merge_AF(annotation_tsv, af_tsv, output_tsv):

    ann_df = pd.read_csv(annotation_tsv, sep="\t", dtype=str)
    af_df  = pd.read_csv(af_tsv, sep="\t", dtype=str)

    # detect ALT1_AF columns automatically
    af_cols = [c for c in af_df.columns if c.endswith("_AF_ALT1")]

    merge_cols = ["CHROM","POS","REF","ALT"]

    af_subset = af_df[merge_cols + af_cols]

    merged = ann_df.merge(
        af_subset,
        on=merge_cols,
        how="left"
    )

    merged.to_csv(output_tsv, sep="\t", index=False)


# ==========================================================
# PIPELINE
# ==========================================================

print("\nStep 1: Extract variants from joint VCF")
csv_to_vcf(INPUT_CSV, JOINT_VCF, FILTERED_VCF)

print("\nStep 2: Annotate variants")
run_snpeff(FILTERED_VCF, ANNOTATED_VCF)

print("\nStep 3: Convert annotated VCF to TSV")
vcf_to_tsv(ANNOTATED_VCF, ANNOTATED_TSV)

print("\nStep 4: Keep first annotation")
keep_first_annotation(ANNOTATED_TSV, FILTERED_ANN)

print("\nStep 5: Merge allele frequencies")
merge_AF(FILTERED_ANN, AF_TSV, FINAL_TSV)

print("\n✔ Pipeline finished")
print("Final file:", FINAL_TSV)
