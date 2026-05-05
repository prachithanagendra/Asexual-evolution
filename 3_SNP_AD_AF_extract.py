#!/usr/bin/env python3

import os
import pysam
import csv

# ==========================
# USER-DEFINED PATHS
# ==========================
JOINT_VCF = "/13_Joint_genotyping/NGF_SK1ref.vcf"
OUTPUT_DIR = "/14_AD_AF_extracted"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================
# SAMPLE GROUP DEFINITIONS
# ==========================
ancestor = "Ancestor2a"

groups = []
for i in range(1, 13):
    group = [f"A{i}", f"A{i+12}", f"A{i+24}", f"A{i+36}"]
    groups.append(group)

# ==========================
# OPEN VCF
# ==========================
vcf = pysam.VariantFile(JOINT_VCF)

# Sanity check
vcf_samples = list(vcf.header.samples)
if ancestor not in vcf_samples:
    raise ValueError("Ancestor2a not found in VCF header")

# ==========================
# HELPER FUNCTION
# ==========================
def extract_sample_metrics(record, sample):
    """
    Returns:
    DP, REF_AD, ALT1_AD, ALT2_AD, REF_AF, ALT1_AF, ALT2_AF
    """
    data = record.samples[sample]

    DP = data.get("DP", 0) or 0

    AD = data.get("AD", None)
    if AD is None:
        return DP, 0, 0, 0, 0.0, 0.0, 0.0

    REF_AD = AD[0] if len(AD) > 0 else 0
    ALT1_AD = AD[1] if len(AD) > 1 else 0
    ALT2_AD = AD[2] if len(AD) > 2 else 0

    total = REF_AD + ALT1_AD + ALT2_AD
    if total > 0:
        REF_AF = REF_AD / total
        ALT1_AF = ALT1_AD / total
        ALT2_AF = ALT2_AD / total
    else:
        REF_AF = ALT1_AF = ALT2_AF = 0.0

    return DP, REF_AD, ALT1_AD, ALT2_AD, REF_AF, ALT1_AF, ALT2_AF

# ==========================
# PROCESS EACH GROUP
# ==========================
for idx, group in enumerate(groups, start=1):

    samples = [ancestor] + group
    samples = [s for s in samples if s in vcf_samples]

    out_tsv = os.path.join(
        OUTPUT_DIR,
        f"Joint_DP_AD_AF_set_{idx}.tsv"
    )

    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")

        # Header
        header = ["CHROM", "POS", "REF", "ALT"]
        for s in samples:
            header.extend([
                f"{s}_DP",
                f"{s}_REF_AD",
                f"{s}_ALT1_AD",
                f"{s}_ALT2_AD",
                f"{s}_REF_AF",
                f"{s}_ALT1_AF",
                f"{s}_ALT2_AF"
            ])
        writer.writerow(header)

        # Iterate over variants
        for record in vcf.fetch():
            alt_str = ",".join(record.alts) if record.alts else "."

            row = [
                record.chrom,
                record.pos,
                record.ref,
                alt_str
            ]

            for s in samples:
                DP, REF_AD, ALT1_AD, ALT2_AD, REF_AF, ALT1_AF, ALT2_AF = \
                    extract_sample_metrics(record, s)

                row.extend([
                    DP,
                    REF_AD,
                    ALT1_AD,
                    ALT2_AD,
                    f"{REF_AF:.4f}",
                    f"{ALT1_AF:.4f}",
                    f"{ALT2_AF:.4f}"
                ])

            writer.writerow(row)

    print(f"Written: {out_tsv}")

print("All files generated successfully.")

