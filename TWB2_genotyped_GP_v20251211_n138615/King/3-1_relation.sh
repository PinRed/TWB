#!/bin/bash
# ==============================================================================
# Pipeline Step: TWB2 Familial Relatedness Check
# Target Dataset: TWB2 cohort
# Inputs:
#   - TWB2_S03_all_QC_clean.{bed,bim,fam} (TWB2 dataset after SNP QC)
# Outputs:
#   - TWB2_S04_qc.{bed,bim,fam} (TWB2 dataset after QC for relatedness check)
#   - TWB2_S04_prune.{bed,bim,fam} (TWB2 dataset after LD pruning)
#   - TWB2_S04_duplicate.kin0 (KING duplicate check results)
#   - TWB2_relatedness_deg2plus.kin0 (KING relatedness check results for 2nd degree relatives)
# Workflow:
#   1. Perform QC on the dataset for relatedness check missingness > 0.01 and  MAF > 0.05
#   2. Perform LD pruning to reduce SNP redundancy 5000kb 1 0.2
#   3. Run King Table filter > 0.0884 to identify related samples (2nd degree relatives and above)
# ==============================================================================

# 從基本QC完的資料中， 取出近親樣本 > 0.0884 (2nd degree relatives)
# 後續利用R做PCAiR
###############################################################
# Step 1: QC
# Step 2: LD pruning
# Step 3: KING duplicate check
# Step 4: KING relatedness check

ARRAY=TWB2

plink2 \
  --bfile ../Primary_SNP_QC/${ARRAY}_S03_all_QC_clean \
  --geno 0.01 \
  --maf 0.05 \
  --make-bed \
  --out ${ARRAY}_S04_qc
  
# LD pruning
plink2 \
  --bfile ${ARRAY}_S04_qc \
  --indep-pairwise 5000kb 1 0.2 \
  --out ${ARRAY}_S04_prune
  
 plink2 \
  --bfile ${ARRAY}_S04_qc \
  --extract ${ARRAY}_S04_prune.prune.in \
  --make-bed \
  --out ${ARRAY}_S04_ldpruned 

# 2. 提取近親與重複樣本 (Kinship > 0.354 表示近乎 Duplicate/MZ twins)
# KING 預設產出 duplicate 報告
king -b ${ARRAY}_S04_ldpruned.bed \
     --duplicate \
     --prefix ${ARRAY}_S04_duplicate
     
# 用King Table filter
plink2 --bfile ${ARRAY}_S04_ldpruned \
  --make-king-table zs cols=id,nsnp,kinship \
  --king-table-filter 0.0884 \
  --out ${ARRAY}_relatedness_deg2plus

# 解壓縮
zstd -d ${ARRAY}_relatedness_deg2plus.kin0.zst \
     -o ${ARRAY}_relatedness_deg2plus.kin0


# Run TWB2_Familial_Relatedness.R