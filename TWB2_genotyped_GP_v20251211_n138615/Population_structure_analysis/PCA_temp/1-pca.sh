twb2fam="../../King/TWB2_S05.fam"
indir="../meta_twb2_pca"
outdir="."

awk '{print $1,$2}' "$twb2fam" > "TWB2_id.txt"


echo "=========================================="
echo "PCA"
echo "=========================================="

plink2 \
  --bfile "${indir}/TWB2_S05_pca" \
  --keep "TWB2_id.txt" \
  --pca 20 approx \
  --out "${outdir}/TWB2_only_pca"


plink2 \
  --bfile "${indir}/TWB2_S05_pca" \
  --remove "TWB2_id.txt" \
  --pca 20 approx \
  --out "${outdir}/TWB2_remove_pca"