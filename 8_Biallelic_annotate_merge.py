#!/usr/bin/env python3

import os
import pandas as pd
import subprocess

# ==========================================================
# USER INPUTS
# ==========================================================

INPUT_CSV  = "/15_BiallelicSNP_analysis/3_wo_fixed_variants/gal6_B_present.tsv"
JOINT_VCF  = "/2_Joint_variantcalling/snps_indels2bp_dp20_biallelic_all_samples.vcf"
AF_TSV     = "/15_BiallelicSNP_analysis/2_DP20_AF_diff_heatmap/02_REF_AF_different/B12_B24.tsv"

WORK_FOLDER  = "/15_BiallelicSNP_analysis/4_Biallelic_annotation"
FINAL_FOLDER = "/15_BiallelicSNP_analysis/5_Biallelic_annotate_merge"

UPDOWN_FOLDER = os.path.join(FINAL_FOLDER, "outside_gene")
OTHER_FOLDER  = os.path.join(FINAL_FOLDER, "inside_gene")

os.makedirs(WORK_FOLDER, exist_ok=True)
os.makedirs(UPDOWN_FOLDER, exist_ok=True)
os.makedirs(OTHER_FOLDER, exist_ok=True)

# ----------------------------------------------------------
# intermediate files
# ----------------------------------------------------------

FILTERED_VCF  = os.path.join(WORK_FOLDER, "gal6_B_filtered.vcf")
ANNOTATED_VCF = os.path.join(WORK_FOLDER, "gal6_B_annotated.vcf")
ANNOTATED_TSV = os.path.join(WORK_FOLDER, "gal6_B_annotated.tsv")

# ----------------------------------------------------------
# final outputs
# ----------------------------------------------------------

UPDOWN_TSV = os.path.join(UPDOWN_FOLDER, "gal6_B_upstream_downstream.tsv")
OTHER_TSV  = os.path.join(OTHER_FOLDER,  "gal6_B_other_variants.tsv")


# ==========================================================
# STEP 1: Extract variants from joint VCF
# ==========================================================

def csv_to_vcf(csv_file, joint_vcf, output_vcf):

    print("Reading variant list")

    df = pd.read_csv(csv_file, sep="\t")

    variant_set = set(
        zip(df["CHROM"].astype(str), df["POS"], df["REF"], df["ALT"])
    )

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
# STEP 2: Run snpEff annotation
# ==========================================================

def run_snpeff(input_vcf, output_vcf):

    print("Running snpEff annotation")

    cmd = f"java -jar snpEff.jar ann SK1_custom {input_vcf} > {output_vcf}"

    subprocess.run(cmd, shell=True, check=True)


# ==========================================================
# STEP 3: Convert annotated VCF → TSV
# ==========================================================

def vcf_to_tsv(input_vcf, output_tsv):

    print("Converting annotated VCF to TSV")

    awk_cmd = r"""awk 'BEGIN {OFS="\t"; print "CHROM","POS","REF","ALT","QUAL","Allele","Consequence","Impact","GeneName","GeneID","FeatureType","FeatureID","TranscriptBioType"} {n = split($6, annots, ","); for(i=1;i<=n;i++){split(annots[i], a, "|"); print $1,$2,$3,$4,$5,a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8]}}'"""

    cmd = f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/ANN\n' {input_vcf} | {awk_cmd} > {output_tsv}"

    subprocess.run(cmd, shell=True, check=True)


# ==========================================================
# STEP 4: Merge AF and apply annotation filtering
# ==========================================================

def merge_and_split(annotation_tsv, af_tsv):

    print("Merging annotation with AF table")

    ann_df = pd.read_csv(annotation_tsv, sep="\t", dtype=str)
    af_df  = pd.read_csv(af_tsv, sep="\t", dtype=str)

    merge_cols = ["CHROM", "POS", "REF", "ALT"]

    af_cols = [c for c in af_df.columns if c.endswith("_AF_ALT1")]

    af_subset = af_df[merge_cols + af_cols]

    merged = ann_df.merge(
        af_subset,
        on=merge_cols,
        how="left"
    )

    outside_groups = []
    inside_groups  = []

    print("Applying annotation rules")

    for _, group in merged.groupby(["CHROM","POS","REF","ALT"]):

        first_cons = str(group.iloc[0]["Consequence"])

        if ("upstream_gene_variant" in first_cons) or ("downstream_gene_variant" in first_cons):

            # keep ALL annotations
            outside_groups.append(group)

        else:

            # keep only first annotation
            inside_groups.append(group.iloc[[0]])

    updown = pd.concat(outside_groups, ignore_index=True) if outside_groups else pd.DataFrame()
    other  = pd.concat(inside_groups, ignore_index=True)  if inside_groups else pd.DataFrame()

    # ------------------------------------------------------
    # Apply ancestor AF filter ONLY to outside_gene variants
    # ------------------------------------------------------

    if "Ancestor2a_AF_ALT1" in updown.columns:

        updown["Ancestor2a_AF_ALT1"] = pd.to_numeric(
            updown["Ancestor2a_AF_ALT1"], errors="coerce"
        )

        updown = updown[
            ~updown["Ancestor2a_AF_ALT1"].between(0.2, 0.8)
        ]

    else:
        print("Warning: Ancestor2a_AF_ALT1 column not found")

    # ------------------------------------------------------
    # Save outputs
    # ------------------------------------------------------

    updown.to_csv(UPDOWN_TSV, sep="\t", index=False)
    other.to_csv(OTHER_TSV, sep="\t", index=False)


# ==========================================================
# PIPELINE
# ==========================================================

print("\nStep 1: Extract variants from joint VCF")
csv_to_vcf(INPUT_CSV, JOINT_VCF, FILTERED_VCF)

print("\nStep 2: Annotate variants")
run_snpeff(FILTERED_VCF, ANNOTATED_VCF)

print("\nStep 3: Convert annotated VCF to TSV")
vcf_to_tsv(ANNOTATED_VCF, ANNOTATED_TSV)

print("\nStep 4: Merge AF and filter annotations")
merge_and_split(ANNOTATED_TSV, AF_TSV)

print("\n✔ Pipeline finished")
print("Outside gene file:", UPDOWN_TSV)
print("Inside gene file:", OTHER_TSV)
