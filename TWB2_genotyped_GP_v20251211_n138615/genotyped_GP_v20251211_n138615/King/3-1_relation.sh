# 從基本QC完的資料中， 取出近親樣本 > 0.0884 (2nd degree relatives)
# 後續利用R做PCAiR
###############################################################
# Step 1: QC
# Step 2: LD pruning
# Step 3: KING duplicate check
# Step 4: KING relatedness check

ARRAY=TWB2

# SNP QC
# geno 0.02
# maf 0.05
plink2 \
  --bfile ../SNP_QC/${ARRAY}_S03_all_QC_clean \
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