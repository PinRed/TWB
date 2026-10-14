#!/bin/bash
# ==============================================================================
# Pipeline Step: Multi-Cohort Data Preparation & Merge & PCA
# Target Datasets: 1KGP, SGDP, and TWB1 cohorts
# 
# Inputs:
#   - 1KGP, SGDP, and TWB1 PLINK binary files (.bed, .bim, .fam)
#
# Outputs:
#   - TWB1_S05_pca.{bed,bim,fam}: Merged, filtered, and LD-pruned dataset
#   - TWB1_S05_pca.eigenvec / eigenval: PCA results (first 20 PCs)
# Do :
#  1. Filter out strand-ambiguous SNPs (A/T and C/G alleles)
#  2. Merge the filtered datasets and catch common SNPs
#  3. Perform LD pruning
#  4. Run PCA
# ==============================================================================

mkdir -p meta_twb1_pca

onekgp="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/1kGP/1kgp_QC_LD"
sdgp="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/SGDP/plink_output/sgdp_hg38_QC"
twb1="../King/TWB1_S05"
outdir="meta_twb1_pca"

# Function to filter out strand-ambiguous SNPs (A/T and C/G alleles) to prevent strand alignment issues during merging
filter_ambiguous() {
  awk '
  !(( ($5=="A") && ($6=="T") ) || \
    ( ($5=="T") && ($6=="A") ) || \
    ( ($5=="C") && ($6=="G") ) || \
    ( ($5=="G") && ($6=="C") )) { print $2 }' "$1"
}

############################################################
# Step 1: Filter Strand-Ambiguous SNPs & Find Intersecting Variants
############################################################

echo "=========================================="
echo "Finding common non-ambiguous SNPs across cohorts"
echo "=========================================="

filter_ambiguous "${onekgp}.bim" | sort -u > 1kg.snps
filter_ambiguous "${sdgp}.bim"   | sort -u > sdgp.snps
filter_ambiguous "${twb1}.bim"   | sort -u > twb1.snps

# Find intersection of SNPs present in all three datasets
comm -12 1kg.snps sdgp.snps > tmp.pos
comm -12 tmp.pos twb1.snps > "${outdir}/common.pos"

echo "Common SNPs count:"
wc -l "${outdir}/common.pos"


############################################################
# Step 2: Extract Common SNPs from Each Dataset
############################################################

echo "=========================================="
echo "Extracting common SNPs per dataset"
echo "=========================================="

plink \
  --bfile "${onekgp}" \
  --extract "${outdir}/common.pos" \
  --make-bed \
  --out "${outdir}/1kg_common"

plink \
  --bfile "${sdgp}" \
  --extract "${outdir}/common.pos" \
  --make-bed \
  --out "${outdir}/sdgp_common"

plink \
  --bfile "${twb1}" \
  --extract "${outdir}/common.pos" \
  --make-bed \
  --out "${outdir}/TWB1_common"


############################################################
# Step 3: Merge Datasets
############################################################

echo "=========================================="
echo "Merging datasets"
echo "=========================================="

# Create input list for PLINK multi-file merge
cat > "${outdir}/merge_common.list" << EOF
${outdir}/sdgp_common.bed ${outdir}/sdgp_common.bim ${outdir}/sdgp_common.fam
${outdir}/TWB1_common.bed ${outdir}/TWB1_common.bim ${outdir}/TWB1_common.fam
EOF

plink \
  --bfile "${outdir}/1kg_common" \
  --merge-list "${outdir}/merge_common.list" \
  --make-bed \
  --out "${outdir}/all_common"


############################################################
# Step 4: LD Pruning
############################################################

echo "=========================================="
echo "LD pruning"
echo "=========================================="

plink \
  --bfile "${outdir}/all_common" \
  --indep-pairwise 1000 1 0.1 \
  --out "${outdir}/TWB1_for_pca"


############################################################
# Step 5: Extract LD-pruned SNPs
############################################################

echo "=========================================="
echo "Extracting LD-pruned SNPs"
echo "=========================================="

plink \
  --bfile "${outdir}/all_common" \
  --extract "${outdir}/TWB1_for_pca.prune.in" \
  --make-bed \
  --out "${outdir}/TWB1_S05_pca"


############################################################
# Step 6: Principal Component Analysis (PCA)
############################################################

echo "=========================================="
echo "Step 6: Running PCA"
echo "=========================================="

plink2 \
  --bfile "${outdir}/TWB1_S05_pca" \
  --pca 20 approx \
  --out "${outdir}/TWB1_S05_pca"


############################################################
# Pipeline Summary & Statistics
############################################################

echo ""
echo "=========================================="
echo "PCA completed successfully!"
echo "=========================================="

echo "Final BED files:"
echo "${outdir}/TWB1_S05_pca.bed"
echo "${outdir}/TWB1_S05_pca.bim"
echo "${outdir}/TWB1_S05_pca.fam"

echo ""
echo "PCA results:"
echo "${outdir}/TWB1_S05_pca.eigenvec"
echo "${outdir}/TWB1_S05_pca.eigenval"

echo "=========================================="
echo "Variant Count Summary"
echo "=========================================="

printf "All common SNPs:     "
wc -l < "${outdir}/common.pos"

printf "LD-pruned SNPs:     "
wc -l < "${outdir}/TWB1_for_pca.prune.in"


############################################################
# Cleanup Intermediate Files
############################################################

echo "=========================================="
echo "Cleanup intermediate files"
echo "=========================================="

# Remove intermediate individual dataset files
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
  "${outdir}/TWB1_common.bed" \
  "${outdir}/TWB1_common.bim" \
  "${outdir}/TWB1_common.fam" \
  "${outdir}/TWB1_common.log" \
  "${outdir}/TWB1_common.nosex" \
  "${outdir}/tmp.pos"

# Remove intermediate merged files
rm -f \
  "${outdir}/merge_common.list" \
  "${outdir}/all_common.bed" \
  "${outdir}/all_common.bim" \
  "${outdir}/all_common.fam" \
  "${outdir}/all_common.log" \
  "${outdir}/all_common.nosex"

# Remove intermediate LD pruning files
rm -f \
  "${outdir}/TWB1_for_pca.prune.in" \
  "${outdir}/TWB1_for_pca.prune.out" \
  "${outdir}/TWB1_for_pca.log" \
  "${outdir}/TWB1_for_pca.nosex"
mkdir -p meta_twb1_pca

onekgp="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/1kGP/1kgp_QC_LD"
sdgp="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/SGDP/plink_output/sgdp_hg38_QC"
twb1="../King/TWB1_S05"
outdir="meta_twb1_pca"

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
filter_ambiguous "${twb1}.bim"   | sort -u > twb1.snps


comm -12 1kg.snps sdgp.snps > tmp.pos
comm -12 tmp.pos twb1.snps > common.pos

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
  --bfile "${twb1}" \
  --extract "${outdir}/common.pos" \
  --make-bed \
  --out "${outdir}/TWB1_common"


echo "=========================================="
echo "Merge datasets"
echo "=========================================="

cat > "${outdir}/merge_common.list" << EOF
${outdir}/sgdp_common.bed ${outdir}/sgdp_common.bim ${outdir}/sgdp_common.fam
${outdir}/TWB1_common.bed ${outdir}/TWB1_common.bim ${outdir}/TWB1_common.fam
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
  --out "${outdir}/TWB1_for_pca"


############################################################
# Step 5. Extract LD-pruned SNPs
############################################################

echo "=========================================="
echo "Extract LD-pruned SNPs"
echo "=========================================="

plink \
  --bfile "${outdir}/all_common" \
  --extract "${outdir}/TWB1_for_pca.prune.in" \
  --make-bed \
  --out "${outdir}/TWB1_S05_pca"


############################################################
# Step 6. PCA
############################################################

echo "=========================================="
echo "Step 6: PCA"
echo "=========================================="

plink2 \
  --bfile "${outdir}/TWB1_S05_pca" \
  --pca 20 approx \
  --out "${outdir}/TWB1_S05_pca"


############################################################
# Finished
############################################################

echo ""
echo "=========================================="
echo "PCA completed!"
echo "=========================================="

echo "Final BED files:"
echo "${outdir}/TWB1_S05_pca.bed"
echo "${outdir}/TWB1_S05_pca.bim"
echo "${outdir}/TWB1_S05_pca.fam"

echo ""
echo "PCA results:"
echo "${outdir}/TWB1_S05_pca.eigenvec"
echo "${outdir}/TWB1_S05_pca.eigenval"

echo "=========================================="
echo "Each Step Number of SNPs"
echo "=========================================="


printf "All common SNPs:     "
wc -l < "${outdir}/common.pos"

printf "LD SNPs:     "
wc -l < "${outdir}/TWB1_for_pca.prune.in"

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
  "${outdir}/TWB1_common.bed" \
  "${outdir}/TWB1_common.bim" \
  "${outdir}/TWB1_common.fam" \
  "${outdir}/TWB1_common.log" \
  "${outdir}/TWB1_common.nosex"\
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
  "${outdir}/TWB1_for_pca.prune.in" \
  "${outdir}/TWB1_for_pca.prune.out" \
  "${outdir}/TWB1_for_pca.log" \
  "${outdir}/TWB1_for_pca.nosex"