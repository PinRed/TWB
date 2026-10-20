#!/usr/bin/env bash
set -euo pipefail

############################################################
# F-SNPQC-01 to F-SNPQC-04
#
# Input:
#   {ARRAY}_S03_all_QC_clean.{bed,bim,fam}
#
# Frozen samples:
#   Population_structure_analysis/{ARRAY}_EAS.id
#
# QC:
#   F-SNPQC-01: Freeze final samples
#   F-SNPQC-02: SNP missingness <= 0.02
#   F-SNPQC-03: Remove monomorphic SNPs (MAF = 0)
#   F-SNPQC-04: HWE P >= 1e-6
#
# Requires PLINK 1.9
############################################################

ARRAY="TWB1"

INPUT="SNP_QC_FILE/TWB1_S06_hg38"
FINAL_SAMPLES="../Population_structure_analysis/${ARRAY}_EAS.id"
OUTDIR="SNP_QC_temp"

mkdir -p "${OUTDIR}"

############################################################
# F-SNPQC-01: Freeze final samples
#
# IMPORTANT:
# --keep normally expects two columns: FID IID.
############################################################

plink \
    --bfile "${INPUT}" \
    --keep "${FINAL_SAMPLES}" \
    --make-bed \
    --out "${OUTDIR}/${ARRAY}_finalsample_base"


############################################################
# F-SNPQC-01.5: Remove duplicate SNP IDs
#
# Keep the first occurrence of each SNP ID.
############################################################


awk '
    seen[$2]++ { print $2 }
' "${OUTDIR}/${ARRAY}_finalsample_base.bim" \
  | sort -u \
  > "${OUTDIR}/${ARRAY}_duplicate_snps.exclude"

printf "Duplicate SNP IDs excluded: %d\n" \
    "$(wc -l < "${OUTDIR}/${ARRAY}_duplicate_snps.exclude")"

plink \
    --bfile "${OUTDIR}/${ARRAY}_finalsample_base" \
    --exclude "${OUTDIR}/${ARRAY}_duplicate_snps.exclude" \
    --make-bed \
    --out "${OUTDIR}/${ARRAY}_deduplicated"

############################################################
# F-SNPQC-02: SNP call rate
############################################################

plink \
    --bfile "${OUTDIR}/${ARRAY}_deduplicated" \
    --geno 0.02 \
    --make-bed \
    --out "${OUTDIR}/${ARRAY}_F_SNPQC02_callrate"


############################################################
# F-SNPQC-03: Calculate allele frequencies
############################################################

plink \
    --bfile "${OUTDIR}/${ARRAY}_F_SNPQC02_callrate" \
    --freq \
    --out "${OUTDIR}/${ARRAY}_F_SNPQC03_freq"

############################################################
# F-SNPQC-03: Exclude monomorphic SNPs (MAF = 0)
#
# PLINK 1.9 .frq columns:
# CHR SNP A1 A2 MAF NCHROBS
# MAF is column 5.
############################################################

awk '
    BEGIN { FS = "[ \t]+" }
    NR > 1 && ($5 + 0) == 0 { print $2 }
' "${OUTDIR}/${ARRAY}_F_SNPQC03_freq.frq" \
  > "${OUTDIR}/${ARRAY}_F_SNPQC03_monomorphic.exclude"

printf "Monomorphic SNPs (MAF = 0) excluded: %d\n" \
    $(wc -l < "${OUTDIR}/${ARRAY}_F_SNPQC03_monomorphic.exclude")

plink \
    --bfile "${OUTDIR}/${ARRAY}_F_SNPQC02_callrate" \
    --exclude "${OUTDIR}/${ARRAY}_F_SNPQC03_monomorphic.exclude" \
    --make-bed \
    --out "${OUTDIR}/${ARRAY}_F_SNPQC03_polymorphic"

############################################################
# F-SNPQC-04: Calculate HWE statistics
############################################################

plink \
    --bfile "${OUTDIR}/${ARRAY}_F_SNPQC03_polymorphic" \
    --hardy \
    --out "${OUTDIR}/${ARRAY}_F_SNPQC04_hwe_stats"

############################################################
# F-SNPQC-04: Exclude SNPs with HWE P < 1e-6
#
# PLINK 1.9 .hwe columns:
# CHR SNP TEST A1 A2 GENO O(HET) E(HET) P
#
# Here TEST=ALL is used because the HWE test is intended
# for the defined sample set as a whole.
############################################################

awk '
    BEGIN { FS = "[ \t]+" }
    NR > 1 && $3 == "ALL" && $9 != "NA" && ($9 + 0) < 1e-6 {
        print $2
    }
' "${OUTDIR}/${ARRAY}_F_SNPQC04_hwe_stats.hwe" \
  | sort -u \
  > "${OUTDIR}/${ARRAY}_F_SNPQC04_hwe.exclude"

############################################################
# F-SNPQC-04: Create final post-HWE fileset
############################################################

plink \
    --bfile "${OUTDIR}/${ARRAY}_F_SNPQC03_polymorphic" \
    --exclude "${OUTDIR}/${ARRAY}_F_SNPQC04_hwe.exclude" \
    --make-bed \
    --out "${OUTDIR}/${ARRAY}_F_SNPQC04_hwe"


############################################################
# Summary
############################################################

echo
echo "F-SNPQC-01 to F-SNPQC-04 completed."
echo "Final files:"
echo "  ${OUTDIR}/${ARRAY}_F_SNPQC04_hwe.bed"
echo "  ${OUTDIR}/${ARRAY}_F_SNPQC04_hwe.bim"
echo "  ${OUTDIR}/${ARRAY}_F_SNPQC04_hwe.fam"
echo
echo "Monomorphic exclusion list:"
echo "  ${OUTDIR}/${ARRAY}_F_SNPQC03_monomorphic.exclude"
echo
echo "HWE exclusion list:"
echo "  ${OUTDIR}/${ARRAY}_F_SNPQC04_hwe.exclude"


rm -rf "${OUTDIR}/${ARRAY}_finalsample_base."* \
       "${OUTDIR}/${ARRAY}_deduplicated."* \
       "${OUTDIR}/${ARRAY}_F_SNPQC02_callrate."* \
       "${OUTDIR}/${ARRAY}_F_SNPQC03_freq."* \
       "${OUTDIR}/${ARRAY}_F_SNPQC03_polymorphic."* \
       "${OUTDIR}/${ARRAY}_F_SNPQC04_hwe_stats."*
