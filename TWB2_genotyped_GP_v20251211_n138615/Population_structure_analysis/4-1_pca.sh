#!/bin/bash
mkdir -p meta_twb2_pca

onekgp="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/1kGP/1kgp_QC_LD"
sdgp="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/SGDP/plink_output/sgdp_hg38_QC"
twb2="../King/TWB2_S05"
outdir="meta_twb2_pca"

filter_ambiguous() {
  awk '
  !(( ($5=="A") && ($6=="T") ) || \
    ( ($5=="T") && ($6=="A") ) || \
    ( ($5=="C") && ($6=="G") ) || \
    ( ($5=="G") && ($6=="C") )) { print $2 }' "$1"
}

# 抓取非 Strand Ambiguous 的位點並排序
filter_ambiguous "${onekgp}.bim" | sort -u > 1kg.snps
filter_ambiguous "${sdgp}.bim"   | sort -u > sdgp.snps
filter_ambiguous "${twb2}.bim"   | sort -u > twb2.snps


comm -12 1kg.snps sdgp.snps > tmp.pos
comm -12 tmp.pos twb2.snps > common.pos

echo "Common SNPs:"
wc -l "${outdir}/common.pos"

plink \
  --bfile "${onekgp}" \
  --extract "${outdir}/common.pos" \
  --make-bed \
  --out "${outdir}/1kg_common"

plink \
  --bfile "${sdgp}" \
  --extract "${outdir}/common.pos" \
  --make-bed \
  --out "${outdir}/sgdp_common"

plink \
  --bfile "${twb2}" \
  --extract "${outdir}/common.pos" \
  --make-bed \
  --out "${outdir}/TWB2_common"


echo "=========================================="
echo "Merge datasets"
echo "=========================================="

cat > "${outdir}/merge_common.list" << EOF
${outdir}/sgdp_common.bed ${outdir}/sgdp_common.bim ${outdir}/sgdp_common.fam
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
echo "LD pruning"
echo "=========================================="

plink \
  --bfile "${outdir}/all_common" \
  --indep-pairwise 1000 1 0.1 \
  --out "${outdir}/TWB2_for_pca"


############################################################
# Step 5. Extract LD-pruned SNPs
############################################################

echo "=========================================="
echo "Extract LD-pruned SNPs"
echo "=========================================="

plink \
  --bfile "${outdir}/all_common" \
  --extract "${outdir}/TWB2_for_pca.prune.in" \
  --make-bed \
  --out "${outdir}/TWB2_S05_pca"


############################################################
# Step 6. PCA
############################################################

echo "=========================================="
echo "Step 6: PCA"
echo "=========================================="

plink2 \
  --bfile "${outdir}/TWB2_S05_pca" \
  --pca 20 approx \
  --out "${outdir}/TWB2_S05_pca"


############################################################
# Finished
############################################################

echo ""
echo "=========================================="
echo "PCA completed!"
echo "=========================================="

echo "Final BED files:"
echo "${outdir}/TWB2_S05_pca.bed"
echo "${outdir}/TWB2_S05_pca.bim"
echo "${outdir}/TWB2_S05_pca.fam"

echo ""
echo "PCA results:"
echo "${outdir}/TWB2_S05_pca.eigenvec"
echo "${outdir}/TWB2_S05_pca.eigenval"

echo "=========================================="
echo "Each Step Number of SNPs"
echo "=========================================="


printf "All common SNPs:     "
wc -l < "${outdir}/common.pos"

printf "LD SNPs:     "
wc -l < "${outdir}/TWB2_for_pca.prune.in"

############################################################
# Cleanup
############################################################

echo "=========================================="
echo "Cleanup intermediate files"
echo "=========================================="

# Step 1: SNP lists
#rm -f \
#  "${outdir}/common.pos"

# Step 2: Common datasets
rm -f \
  "${outdir}/1kg_common.bed" \
  "${outdir}/1kg_common.bim" \
  "${outdir}/1kg_common.fam" \
  "${outdir}/1kg_common.log" \
  "${outdir}/1kg_common.nosex" \
  "${outdir}/sdgp_common.bed" \
  "${outdir}/sdgp_common.bim" \
  "${outdir}/sdgp_common.fam" \
  "${outdir}/sdgp_common.log" \
  "${outdir}/sdgp_common.nosex" \
  "${outdir}/TWB2_common.bed" \
  "${outdir}/TWB2_common.bim" \
  "${outdir}/TWB2_common.fam" \
  "${outdir}/TWB2_common.log" \
  "${outdir}/TWB2_common.nosex"\
  "${outdir}/tmp.pos" \

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