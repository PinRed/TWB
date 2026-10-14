ARRAY=TWB2
twb2_prelim="../King/${ARRAY}_S04_batch_SNPs_removed"
eas_id="${ARRAY}_EAS.id"

#########################################
# 抓EAS族群的個體，
# 先做SNP QC、LD pruning，再做KING，抓出親緣關係
# input : TWB2 做完 sample QC 、 Preliminary SNP QC 與 batch effect SNPs removal 後的資料
#########################################

echo "=========================================="
echo "SNP QC and LD pruning"
echo "=========================================="

plink2 --bfile ${twb2_prelim} \
  --keep ${eas_id} \
  --autosome \
  --geno 0.01 \
  --maf 0.05 \
  --indep-pairwise 5000kb 0.3 \
  --out ${ARRAY}_EAS_pcair_prune


plink2 --bfile ${twb2_prelim} \
  --keep ${eas_id} \
  --extract ${ARRAY}_EAS_pcair_prune.prune.in \
  --make-king-table zs cols=id,nsnp,kinship \
  --king-table-filter 0.02209709 \
  --out ${ARRAY}_EAS_pcair_kin

plink2 --zst-decompress \
  ${ARRAY}_EAS_pcair_kin.kin0.zst \
  ${ARRAY}_EAS_pcair_kin.kin0
