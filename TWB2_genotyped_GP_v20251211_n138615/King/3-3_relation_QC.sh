#!/bin/bash
# ==============================================================================
# Pipeline Step: TWB2 Relatedness Check
# Target Dataset: TWB2 cohort
# Inputs:
#   - TWB2_S03_all_QC_clean.{bed,bim,fam} (TWB2 dataset after SNP QC)
# Outputs:
#   - TWB2_S04_batch_SNPs_removed.{bed,bim,fam} (TWB2 dataset after removing batch level SNPs)
#   - TWB2_S05.{bed,bim,fam} (TWB2 dataset after removing related samples)
# Workflow:
#   1. Remove batch level SNPs and related samples from the dataset
#   2. Create a new dataset after removing the failed samples
# ==============================================================================

# 過濾Batch SNP missingness、allene-frequency 以及 related samples
ARRAY=TWB2

tail -n +2 ${ARRAY}_filter.relatedsamples | \
awk '{print $1, $1}' \
> ${ARRAY}_filter.remove

# Catch common SNPs and remove related samples
plink2 \
  --bfile ../Primary_SNP_QC/${ARRAY}_S03_all_QC_clean \
  --extract ../Batch-level_SNP_QC/intersect_snps.id \
  --make-bed \
  --out ${ARRAY}_S04_batch_SNPs_removed

plink2 \
  --bfile ${ARRAY}_S04_batch_SNPs_removed \
  --remove ${ARRAY}_filter.remove \
  --set-all-var-ids '@:#:$r:$a' \
  --make-bed \
  --out ${ARRAY}_S05

  rm -rf ${ARRAY}_S04_qc.* \
       ${ARRAY}_S04_prune.* \
       ${ARRAY}_S04_ldpruned.* \
       ${ARRAY}_S04_duplicate.* \
       ${ARRAY}_relatedness_deg2plus.kin0.zst