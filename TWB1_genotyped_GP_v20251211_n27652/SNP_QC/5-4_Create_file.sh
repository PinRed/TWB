#!/usr/bin/env bash
mkdir -p SNP_QC_FILE

ARRAY=TWB1
Input_BFILE="SNP_QC_temp/${ARRAY}_F_SNPQC04_hwe"
Exclude_SNPs="TWB1_F_SNPQC05_batch_effect/problematic_variants_to_exclude.txt"
OUTDIR="SNP_QC_FILE/"

plink2 \
  --bfile "${Input_BFILE}" \
  --exclude "${Exclude_SNPs}" \
  --set-all-var-ids '@:#:$r:$a' \
  --make-bed \
  --out "${OUTDIR}${ARRAY}_F_SNPQC05_clean"

rm -rf SNP_QC_temp/

mkdir -p Freq
OUT="Freq/"

META="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/metadata/igsr_samples.tsv"
BFILE="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/1kgp/1kgp_hg19_QC_LD"

# 1. 擷取 EAS sample ID
awk -F '\t' 'NR>1 && $6=="EAS" {print $1}' \
  "$META" > "${OUT}EAS_ID.txt"

# 2. 與 PLINK .fam IID 比對，輸出真正的 FID IID
awk '
NR==FNR {eas[$1]=1; next}
$2 in eas {print $1, $2}
' "${OUT}EAS_ID.txt" "${BFILE}.fam" \
  > "${OUT}Onekg_EAS_sample_list.txt"

# 3. 確認樣本數
wc -l "${OUT}Onekg_EAS_sample_list.txt"

# 4. 計算 EAS allele frequency
plink \
  --bfile "$BFILE" \
  --keep "${OUT}Onekg_EAS_sample_list.txt" \
  --freq \
  --out "${OUT}1KGP_EAS_freq"

plink \
  --bfile "${OUTDIR}${ARRAY}_F_SNPQC05_clean" \
  --freq \
  --out "${OUT}TWB_freq"