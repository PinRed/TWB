#!/usr/bin/env bash
mkdir -p SNP_QC_FILE

ARRAY=TWB1
Input_BFILE="SNP_QC_temp/${ARRAY}_F_SNPQC04_hwe"
Exclude_SNPs="TWB1_F_SNPQC05_batch_effect/problematic_variants_to_exclude.txt"
OUTDIR="SNP_QC_FILE/"

plink2 \
  --bfile "${Input_BFILE}" \
  --exclude "${Exclude_SNPs}" \
  --make-bed \
  --out "${OUTDIR}${ARRAY}_F_SNPQC05_clean"

plink2 \
    --bfile "${OUTDIR}${ARRAY}_F_SNPQC05_clean" \
    --make-pgen \
    --out "${OUTDIR}${ARRAY}_F_SNPQC06_base"

rm -rf SNP_QC_temp/