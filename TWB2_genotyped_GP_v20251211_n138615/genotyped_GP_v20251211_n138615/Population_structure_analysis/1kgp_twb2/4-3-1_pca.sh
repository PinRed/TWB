#!/bin/bash
# Only common SNPs between 1KG and TWB2 are used for PCA analysis.

set -e

############################################################
# Input
############################################################

# onekg_all="../../../1kGP_SGDP_prep/1kGP/1kgp_QC_LD"
# sgdp_all="../../../1kGP_SGDP_prep/SGDP/plink_output/sgdp_QC"
# twb2_all="../King/TWB2_S05"

onekg_all="harmonized_common/1kgp_common"
twb2_all="harmonized_common/twb2_common"

outdir="."

mkdir -p "${outdir}"


############################################################
# Step 1. Find intersection SNPs
# 1KG ∩ SGDP ∩ TWB2
############################################################

echo "=========================================="
echo "Step 1: Find common SNPs"
echo "=========================================="

awk '{print $2}' "${onekg_all}.bim" | sort -u \
  > "${outdir}/1kg.snps"

awk '{print $2}' "${twb2_all}.bim" | sort -u \
  > "${outdir}/twb2.snps"

comm -12 \
  "${outdir}/1kg.snps" \
  "${outdir}/twb2.snps" \
  > "${outdir}/1kg_TWB2_common.id"


echo "Common SNPs:"
wc -l "${outdir}/1kg_TWB2_common.id"


############################################################
# Step 2. Extract common SNPs
############################################################

echo "=========================================="
echo "Step 2: Extract common SNPs"
echo "=========================================="

plink \
  --bfile "${onekg_all}" \
  --extract "${outdir}/1kg_TWB2_common.id" \
  --make-bed \
  --out "${outdir}/1kg_common"


plink \
  --bfile "${twb2_all}" \
  --extract "${outdir}/1kg_TWB2_common.id" \
  --make-bed \
  --out "${outdir}/TWB2_common"


############################################################
# Step 3. Merge 1KG + SGDP + TWB2
############################################################

echo "=========================================="
echo "Step 3: Merge datasets"
echo "=========================================="

cat > "${outdir}/merge_common.list" << EOF
${outdir}/1kg_common.bed ${outdir}/1kg_common.bim ${outdir}/1kg_common.fam
${outdir}/TWB2_common.bed ${outdir}/TWB2_common.bim ${outdir}/TWB2_common.fam
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
echo "Step 4: LD pruning"
echo "=========================================="

plink \
  --bfile "${outdir}/all_common" \
  --indep-pairwise 1000 1 0.1 \
  --out "${outdir}/TWB2_for_pca"


############################################################
# Step 5. Extract LD-pruned SNPs
############################################################

echo "=========================================="
echo "Step 5: Extract LD-pruned SNPs"
echo "=========================================="

plink \
  --bfile "${outdir}/all_common" \
  --extract "${outdir}/TWB2_for_pca.prune.in" \
  --make-bed \
  --out "${outdir}/TWB2_pca"


############################################################
# Step 6. PCA
############################################################

echo "=========================================="
echo "Step 6: PCA"
echo "=========================================="

plink2 \
  --bfile "${outdir}/TWB2_pca" \
  --pca 20 approx \
  --out "${outdir}/TWB2_pca"


############################################################
# Finished
############################################################

echo ""
echo "=========================================="
echo "PCA completed!"
echo "=========================================="

echo "Final BED files:"
echo "${outdir}/TWB2_pca.bed"
echo "${outdir}/TWB2_pca.bim"
echo "${outdir}/TWB2_pca.fam"

echo ""
echo "PCA results:"
echo "${outdir}/TWB2_pca.eigenvec"
echo "${outdir}/TWB2_pca.eigenval"

echo "=========================================="
echo "Each Step Number of SNPs"
echo "=========================================="

printf "1KG SNPs:              "
wc -l < "${outdir}/1kg.snps"

printf "TWB1 SNPs:             "
wc -l < "${outdir}/twb2.snps"

printf "All common SNPs:     "
wc -l < "${outdir}/1kg_TWB2_common.id"

printf "LD SNPs:     "
wc -l < "${outdir}/TWB2_for_pca.prune.in"

############################################################
# Cleanup
############################################################

echo "=========================================="
echo "Cleanup intermediate files"
echo "=========================================="

# Step 1: SNP lists
rm -f \
  "${outdir}/1kg.snps" \
  "${outdir}/twb2.snps" \
  "${outdir}/1kg_sgdp_common.id" \
  "${outdir}/common.id"

# Step 2: Common datasets
rm -f \
  "${outdir}/1kg_common.bed" \
  "${outdir}/1kg_common.bim" \
  "${outdir}/1kg_common.fam" \
  "${outdir}/1kg_common.log" \
  "${outdir}/1kg_common.nosex" \
  "${outdir}/TWB2_common.bed" \
  "${outdir}/TWB2_common.bim" \
  "${outdir}/TWB2_common.fam" \
  "${outdir}/TWB2_common.log" \
  "${outdir}/TWB2_common.nosex"

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
  "${outdir}/TWB2_for_pca.prune.in" \
  "${outdir}/TWB2_for_pca.prune.out" \
  "${outdir}/TWB2_for_pca.log" \
  "${outdir}/TWB2_for_pca.nosex"

# Step 5/6: PCA input genotype files
#  rm -f \
#   "${outdir}/TWB2_pca.bed" \
#   "${outdir}/TWB2_pca.bim" \
#   "${outdir}/TWB2_pca.fam" \
#   "${outdir}/TWB2_pca.log" \
#   "${outdir}/TWB2_pca.nosex"
