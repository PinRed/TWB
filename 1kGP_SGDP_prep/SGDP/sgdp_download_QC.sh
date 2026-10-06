# SGDP (Simons Genome Diversity Project)
cd ../SGDP
# 下載 SGDP 樣本 Metadata
wget https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/vcf_variants/vcfs.variants.public_samples.279samples.file.list 

# 下載 SGDP VCF 檔案庫 (EBI / Harvard Reich Lab Archive)
wget https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/vcf_variants/vcfs.variants.public_samples.279samples.tar
# 註：也可從 Harvard Reich Lab 官方空間獲取 PLINK/VCF 格式

# 提取 Region 為 Oceania, EastAsia, SouthAsia 的 Sample ID
awk -F'\t' '$7~/(Oceania|EastAsia|SouthAsia)/ {print $1}' SGDP_metadata.279public.21signedLetter.44Fan.samples.txt > SGDP_SAMPLE_LIST.txt

# 計算 .vcf.gz 主檔案數量
ls *.variants.vcf.gz | wc -l

# =================================================
# 刪除舊的 4KB 空檔案
rm -rf sliced_vcfs/
mkdir -p sliced_vcfs

for vcf in data/*.annotated.nh2.variants.vcf.gz; do
  base_vcf=$(basename "$vcf")
  sample_in_vcf=$(bcftools query -l "$vcf")
  
  if grep -qx "$sample_in_vcf" SGDP_SAMPLE_LIST.txt; then
    echo "Processing matched sample: ${sample_in_vcf} (${base_vcf})..."
    
    # chr 的 1-22 染色體名稱
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
                  -O z -o "sliced_vcfs/sliced_${base_vcf}" "$vcf" --threads 4 && \
    bcftools index -t "sliced_vcfs/sliced_${base_vcf}"
  fi
done

# ============================================================
# 合併 VCF
VCF_LIST="merged_output/vcf_list.txt"
ls sliced_vcfs/sliced_*.vcf.gz > "$VCF_LIST"

echo "=================================================="
echo "Step 1: Merging sliced VCFs into Multi-sample VCF"
echo "=================================================="
# 透過 merge 將所有單樣本合為一檔，bcftools 會自動補齊 0/0 (Hom-Ref)
bcftools merge -l "$VCF_LIST" \
               -O z -o merged_output/sgdp_merged.vcf.gz \
               --threads 4

bcftools index -t merged_output/sgdp_merged.vcf.gz

# ==============================================================
# 轉成 PLINK 格式

plink2 \
    --vcf merged_output/sgdp_merged.vcf.gz \
    --autosome \
    --snps-only just-acgt \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --make-bed \
    --out plink_output/sgdp_autosome_biallelic_snps
    
# Filter MAF < 0.01
plink2 \
  --bfile plink_output/sgdp_autosome_biallelic_snps \
  --maf 0.01 \
  --make-bed \
  --out plink_output/sgdp_QC