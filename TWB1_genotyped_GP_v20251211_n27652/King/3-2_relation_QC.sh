# 過濾Batch SNP missingness、allene-frequency 以及 related samples
ARRAY=TWB1

tail -n +2 ${ARRAY}_filter.relatedsamples | \
awk '{print $1, $1}' \
> ${ARRAY}_filter.remove

# Catch common SNPs and remove related samples
plink2 \
  --bfile ../Primary_SNP_QC/${ARRAY}_S03_all_QC_clean \
  --extract ../SNP_QC_Batch/intersect_snps.id \
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