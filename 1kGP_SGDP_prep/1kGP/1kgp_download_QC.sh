# cd metadata
# # 下載 IGSR 樣本清單 (igsr_samples.tsv)
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.sequence.index -O igsr_samples.tsv 2>/dev/null || \
# curl -o igsr_samples.tsv https://www.internationalgenome.org/api/beta/sample/download

# # 下載 3,202 樣本家系資訊 (1kGP.3202_samples.pedigree_info.txt )
# wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20201028_2504_genomes_3202_samples_pedigree_info.txt -O 1kGP.3202_samples.pedigree_info.txt
# cd ..

cd 1kGP
# 批次下載 Chr1 到 Chr22 的 VCF 及其索引檔
for chr in {1..22}; do
    echo "=== Processing chr${chr} ==="
    
    BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
    VCF_FILE="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chr}.filtered.shapeit2-duohmm-phased.vcf.gz"
    
    # 1. Download VCF and Index
    wget -c --timeout=30 --tries=0 --retry-connrefused "${BASE_URL}/${VCF_FILE}"
    wget -c --timeout=30 --tries=0 --retry-connrefused "${BASE_URL}/${VCF_FILE}.tbi"
    
    # 2. Filter for unrelated samples and biallelic SNPs
    bcftools view \
        -m2 -M2 \
        -v snps \
        -S ../metadata/keep_unrelated_samples.txt \
        -Oz \
        -o "chr${chr}.unrelated.biallelic-snps.vcf.gz" \
        "${VCF_FILE}"
        
    # 3. Index the newly filtered VCF
    bcftools index -t "chr${chr}.unrelated.biallelic-snps.vcf.gz"
    
    # Optional: Delete original large raw file to save disk space
    # rm "${VCF_FILE}" "${VCF_FILE}.tbi"
done

# ==============================================================
# 轉成 PLINK 格式
for chr in {1..22}; do
    echo "=== Converting chr${chr} to PLINK2 format ==="
    plink2 \
        --vcf chr${chr}.unrelated.biallelic-snps.vcf.gz \
        --make-pgen \
        --out chr${chr}.unrelated.biallelic-snps
done
# ===========================================================
# 合併 Chr1-22 PLINK 檔案
for chr in {2..22}; do
    echo "chr${chr}.unrelated.biallelic-snps" >> merge_list.txt
done

# 使用 PLINK2 合併全部染色體
plink2 \
    --pfile chr1.unrelated.biallelic-snps \
    --pmerge-list merge_list.txt \
    --make-bed \
    --out 1kgp_all_chrs.unrelated.biallelic-snps

# 轉換HighLD格式
sed 's/^chr//' \
    HIGH_LD_REGION_BED_hg38.txt \
    > 1kgp_HIGH_LD_REGION_BED_hg38.txt

#===========================================================
# 2. 準備高 LD 區域區間 (high_ld_regions.txt，格式：CHR POS_START POS_END ID)
# 3. 剔除高 LD 區域並提取無親緣樣本
plink2 \
    --bfile 1kgp_all_chrs.unrelated.biallelic-snps \
    --keep ../metadata/keep_unrelated_samples.txt \
    --exclude bed0 ../metadata/HIGH_LD_REGION_BED_all_builds/1kgp_HIGH_LD_REGION_BED_hg38.txt
    --make-bed \
    --out 1kgp_QC

# Filter MAF < 0.01
plink2 \
  --bfile 1kgp_QC \
  --maf 0.01 \
  --make-bed \
  --out 1kgp_QC_LD