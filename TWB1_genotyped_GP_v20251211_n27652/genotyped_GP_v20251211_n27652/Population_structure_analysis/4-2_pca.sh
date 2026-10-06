#!/bin/bash

onekgp="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/1kGP/1kgp_QC_LD"
twb1="./TWB1_S05_hg38"
outdir="."

cut -f2 "${onekgp}.bim" | sort -u > 1kg.snps
cut -f2 "${twb1}.bim" | sort -u > twb1.snps

comm -12 1kg.snps twb1.snps > common.pos



plink \
  --bfile "${onekgp}" \
  --extract "./common.pos" \
  --make-bed \
  --out "${outdir}/1kg_common"


plink \
  --bfile "${twb1}" \
  --extract "${outdir}/common.pos" \
  --make-bed \
  --out "${outdir}/twb1_common"


echo "=========================================="
echo "Merge datasets"
echo "=========================================="

cat > "${outdir}/merge_common.list" << EOF
${outdir}/twb1_common.bed ${outdir}/twb1_common.bim ${outdir}/twb1_common.fam
EOF



plink \
  --bfile "${outdir}/1kg_common" \
  --merge-list "${outdir}/merge_common.list" \
  --make-bed \
  --out "${outdir}/all_common"


############################################################
# Step 4. LD pruning
############################################################

echo "=========================================="
echo "LD pruning"
echo "=========================================="

plink \
  --bfile "${outdir}/all_common" \
  --indep-pairwise 1000 1 0.1 \
  --out "${outdir}/twb1_for_pca"


############################################################
# Step 5. Extract LD-pruned SNPs
############################################################

echo "=========================================="
echo "Extract LD-pruned SNPs"
echo "=========================================="

plink \
  --bfile "${outdir}/all_common" \
  --extract "${outdir}/twb1_for_pca.prune.in" \
  --make-bed \
  --out "${outdir}/twb1_pca"


############################################################
# Step 6. PCA
############################################################

echo "=========================================="
echo "Step 6: PCA"
echo "=========================================="

plink2 \
  --bfile "${outdir}/twb1_pca" \
  --pca 20 approx \
  --out "${outdir}/twb1_pca"


############################################################
# Finished
############################################################

echo ""
echo "=========================================="
echo "PCA completed!"
echo "=========================================="

echo "Final BED files:"
echo "${outdir}/twb1_pca.bed"
echo "${outdir}/twb1_pca.bim"
echo "${outdir}/twb1_pca.fam"

echo ""
echo "PCA results:"
echo "${outdir}/twb1_pca.eigenvec"
echo "${outdir}/twb1_pca.eigenval"

echo "=========================================="
echo "Each Step Number of SNPs"
echo "=========================================="


printf "All common SNPs:     "
wc -l < "${outdir}/common.pos"

printf "LD SNPs:     "
wc -l < "${outdir}/twb1_for_pca.prune.in"

############################################################
# Cleanup
############################################################

echo "=========================================="
echo "Cleanup intermediate files"
echo "=========================================="

# Step 1: SNP lists
rm -f \
  "${outdir}/common.pos"

# Step 2: Common datasets
rm -f \
  "${outdir}/1kg_common.bed" \
  "${outdir}/1kg_common.bim" \
  "${outdir}/1kg_common.fam" \
  "${outdir}/1kg_common.log" \
  "${outdir}/1kg_common.nosex" \
  "${outdir}/twb1_common.bed" \
  "${outdir}/twb1_common.bim" \
  "${outdir}/twb1_common.fam" \
  "${outdir}/twb1_common.log" \
  "${outdir}/twb1_common.nosex"

# Step 3: Merge files
rm -f \
  "${outdir}/merge_common.list" \
  "${outdir}/all_common.bed" \
  "${outdir}/all_common.bim" \
  "${outdir}/all_common.fam" \
  "${outdir}/all_common.log" \
  "${outdir}/all_common.nosex"

# Step 4: LD pruning files
rm -f \
  "${outdir}/twb1_for_pca.prune.in" \
  "${outdir}/twb1_for_pca.prune.out" \
  "${outdir}/twb1_for_pca.log" \
  "${outdir}/twb1_for_pca.nosex"