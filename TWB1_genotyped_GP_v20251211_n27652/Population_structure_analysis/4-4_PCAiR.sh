#!/bin/bash
# ==============================================================================
# Pipeline Step: EAS Sample Extraction, SNP QC, LD Pruning & Kinship Analysis
# Target Dataset: TWB2 cohort (EAS)
# 
# Inputs:
#   - TWB1_S04_batch_SNPs_removed.{bed,bim,fam} 
#     (TWB1 post-sample QC, preliminary SNP QC, and batch effect SNP removal)
#   - TWB1_EAS.id (List of target East Asian individuals to keep)
#
# Outputs:
#   - TWB1_EAS_pcair_prune.{prune.in,prune.out}: LD-pruned SNP list for PC-AIR
#   - TWB1_EAS_pcair_kin.kin0.zst: Compressed KING kinship coefficient table
#   - TWB1_EAS_pcair_kin.kin0: Decompressed KING kinship table filtered at > 3rd-degree ($\ge 0.02209709$)
#
# Do:
#   1. Extract EAS individuals, filter autosomes, apply missingness (geno 0.01) and MAF (0.05) filters, and perform LD pruning 5000kb 0.3
#   2. Compute genetic relationship coefficients (KING) using the pruned SNP set
#   3. Filter related pairs using a 3rd-degree relationship threshold ($0.02209709$)
#   4. Decompress the resulting KING output file
# ==============================================================================

ARRAY=TWB1
twb1_prelim="../King/${ARRAY}_S04_batch_SNPs_removed"
eas_id="${ARRAY}_EAS.id"

#########################################
# 抓EAS族群的個體，
# 先做SNP QC、LD pruning，再做KING，抓出親緣關係
# input : TWB1 做完 sample QC 、 Preliminary SNP QC 與 batch effect SNPs removal 後的資料
#########################################

echo "=========================================="
echo "SNP QC and LD pruning"
echo "=========================================="

plink2 --bfile ${twb1_prelim} \
  --keep ${eas_id} \
  --autosome \
  --geno 0.01 \
  --maf 0.05 \
  --indep-pairwise 5000kb 0.3 \
  --out ${ARRAY}_EAS_pcair_prune


plink2 --bfile ${twb1_prelim} \
  --keep ${eas_id} \
  --extract ${ARRAY}_EAS_pcair_prune.prune.in \
  --make-king-table zs cols=id,nsnp,kinship \
  --king-table-filter 0.02209709 \
  --out ${ARRAY}_EAS_pcair_kin

plink2 --zst-decompress \
  ${ARRAY}_EAS_pcair_kin.kin0.zst \
  ${ARRAY}_EAS_pcair_kin.kin0
