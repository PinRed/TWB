#!/bin/bash
set -Eeuo pipefail

mkdir -p 1kgp_hg38
cd 1kgp_hg38

MAX_ATTEMPTS=10
BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/"

get_vcf_name() {
    echo "CCDG_14151_B01_GRM_WGS_2020-08-05_chr${1}.filtered.shapeit2-duohmm-phased.vcf.gz"
}

log() {
    echo "[$(date '+%F %T')] $*"
}

download_and_validate() {
    local chr="$1"
    local vcf
    local idx
    local attempt

    vcf="$(get_vcf_name "$chr")"
    idx="${vcf}.tbi"

    for attempt in $(seq 1 "$MAX_ATTEMPTS"); do
        log "chr${chr}: attempt ${attempt}/${MAX_ATTEMPTS}"

        # 直接下載正式 VCF 檔案，不使用 .part
        rm -f "$vcf"

        if ! wget --timeout=60 --tries=3 --waitretry=5 \
            -O "$vcf" "${BASE_URL}${vcf}"; then
            log "WARNING: chr${chr} VCF download failed"
            rm -f "$vcf"
            continue
        fi

        # 完整讀取 VCF，檢查資料區塊是否損毀
        if ! bcftools view -Ou "$vcf" >/dev/null; then
            log "WARNING: chr${chr} VCF full-read check failed"
            rm -f "$vcf"
            continue
        fi

        # 下載索引，但不檢查索引內容
        wget --timeout=60 --tries=3 --waitretry=5 \
            -O "$idx" "${BASE_URL}${idx}" || \
            log "WARNING: chr${chr} index download failed"

        log "chr${chr}: VCF validation passed"
        return 0
    done

    log "ERROR: chr${chr} failed after ${MAX_ATTEMPTS} attempts"
    return 1
}


# 批次下載 Chr1 到 Chr22 的 VCF 及其索引檔
for chr in {1..22}; do
    echo "=== Processing chr${chr} ==="
    vcf="$(get_vcf_name "$chr")"

    # 1. Download VCF and Index
    download_and_validate "$chr"
    
    # 2. Filter for unrelated samples and biallelic SNPs
    bcftools view \
        -m2 -M2 \
        -v snps \
        -S ../../metadata/keep_unrelated_samples.txt \
        -Oz \
        -o "chr${chr}.unrelated.biallelic-snps.vcf.gz" \
        "${vcf}"
        
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


#===========================================================
# 2. 準備高 LD 區域區間 (high_ld_regions.txt，格式：CHR POS_START POS_END ID)
# 3. 剔除高 LD 區域並提取無親緣樣本
plink2 \
    --bfile 1kgp_all_chrs.unrelated.biallelic-snps \
    --keep ../../metadata/keep_unrelated_samples.txt \
    --exclude bed0 ../../metadata/HIGH_LD_REGION_BED_all_builds/1kgp_HIGH_LD_REGION_BED_hg38.txt \
    --make-bed \
    --out 1kgp_QC

# Filter MAF < 0.01
plink2 \
  --bfile 1kgp_QC \
  --maf 0.01 \
  --make-bed \
  --out 1kgp_QC_LD