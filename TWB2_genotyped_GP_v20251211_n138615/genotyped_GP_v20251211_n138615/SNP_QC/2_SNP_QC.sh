ARRAY=TWB2
# Run TWB-SNP_QC.R
# {ARRAY}_PSNPQC01_structural_keep.snplist 


plink2 --bfile ../Sample_QC/${ARRAY}_S03_QC \
  --extract ${ARRAY}_PSNPQC01_structural_keep.snplist \
  --make-bed \
  --out ${ARRAY}_S03_QC_PSNPQC01_structural

plink2 --bfile ${ARRAY}_S03_QC_PSNPQC01_structural --freq --out "${ARRAY}_PSNPQC02_freq"

awk 'NR > 1 && $5 == 0 {print $2}' \
  ${ARRAY}_PSNPQC02_freq.afreq \
  > "${ARRAY}_PSNPQC02_monomorphic.snplist"
  
# 由 frequency output 僅建立 MAF == 0 的 SNP list
plink2 --bfile ${ARRAY}_PSNPQC01_structural \
  --exclude ${ARRAY}_PSNPQC02_monomorphic.snplist \
  --make-bed \
  --out "${ARRAY}_PSNPQC02_polymorphic"


# F_MISS > 0.02 SNPs remove
plink2 --bfile ${ARRAY}_PSNPQC02_polymorphic \
  --geno 0.02 \
  --make-bed \
  --out ${ARRAY}_S03_all_QC_clean