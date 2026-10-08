#!/bin/bash

set -euo pipefail

mkdir -p chr_pos

onekgp="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/1kGP/1kgp_QC_LD"
sgdp="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/SGDP/plink_output/sgdp_hg38_QC"
twb2="../../King/TWB2_S05"

plink2 \
  --bfile "${onekgp}" \
  --set-all-var-ids '@:#' \
  --make-bed \
  --out chr_pos/1kgp_pos

plink2 \
  --bfile "${sgdp}" \
  --set-all-var-ids '@:#' \
  --make-bed \
  --out chr_pos/sgdp_pos

plink2 \
  --bfile "${twb2}" \
  --set-all-var-ids '@:#' \
  --make-bed \
  --out chr_pos/twb2_pos

cd chr_pos

# 因為1kgp、sgdp_hg38有duplicated SNPs, 所以要先用去除
cut -f2 1kgp_pos.bim |
sort |
uniq -d > 1kgp.duplicate.pos

cut -f2 sgdp_pos.bim |
sort |
uniq -d > sgdp.duplicate.pos

plink2 \
    --bfile 1kgp_pos \
    --exclude 1kgp.duplicate.pos \
    --make-bed \
    --out 1kgp_unique

plink2 \
    --bfile sgdp_pos \
    --exclude sgdp.duplicate.pos \
    --make-bed \
    --out sgdp_unique

# rm 1kgp.duplicate.pos
rm -rf 1kgp_pos.*
for f in ./1kgp_unique.*; do
    mv "$f" "./1kgp_pos.${f#./1kgp_unique.}"
done

# rm sgdp.duplicate.pos
rm -rf sgdp_pos.*
for f in ./sgdp_unique.*; do
    mv "$f" "./sgdp_pos.${f#./sgdp_unique.}"
done


cut -f2 1kgp_pos.bim | sort -u > 1kgp.pos
cut -f2 sgdp_pos.bim  | sort -u > sgdp.pos
cut -f2 twb2_pos.bim  | sort -u > twb2.pos

comm -12 1kgp.pos sgdp.pos > tmp.pos
comm -12 tmp.pos twb2.pos > common.pos

# rm tmp.pos

echo "Common CHR:POS:"
wc -l common.pos