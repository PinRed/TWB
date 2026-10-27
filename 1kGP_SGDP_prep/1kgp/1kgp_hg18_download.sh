#!/bin/bash
set -Eeuo pipefail

BASE_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
OUTDIR="phase3_1kgp_hg19"
MAX_ATTEMPTS=10

mkdir -p "$OUTDIR"
cd "$OUTDIR"

log() {
    echo "[$(date '+%F %T')] $*"
}

get_vcf_name() {
    echo "ALL.chr${1}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
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

# Step 1: Download all chromosomes and validate each VCF
for chr in {1..22}; do
    download_and_validate "$chr"
done

# Step 2: Final verification of all VCF files only
log "Final VCF verification"

for chr in {1..22}; do
    vcf="$(get_vcf_name "$chr")"

    [[ -s "$vcf" ]] || {
        log "ERROR: Missing or empty VCF: $vcf"
        exit 1
    }

    gzip -t "$vcf"
    bcftools view -Ou "$vcf" >/dev/null

    log "chr${chr}: VCF passed final verification"
done


log "All 22 VCF files passed validation"

# Step 3: Rebuild indexes when needed by downstream processing
for chr in {1..22}; do
    vcf="$(get_vcf_name "$chr")"
    bcftools view \
        -m2 -M2 \
        -v snps \
        -Oz \
        -o "chr${chr}.unrelated.biallelic-snps.vcf.gz" \
        "${vcf}"
    bcftools index -f -t "$vcf"
done

log "Download and VCF validation completed successfully"

for chr in {1..22}; do
    echo "=== Converting chr${chr} to PLINK2 format ==="

    plink2 \
        --vcf chr${chr}.unrelated.biallelic-snps.vcf.gz \
        --set-missing-var-ids '@:#:$r:$a' \
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
    --exclude bed0 /SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/metadata/HIGH_LD_REGION_BED_all_builds/HIGH_LD_REGION_BED_hg37.txt \
    --make-bed \
    --out 1kgp_hg19_QC

# Filter MAF < 0.01
plink2 \
  --bfile 1kgp_hg19_QC \
  --maf 0.01 \
  --make-bed \
  --out 1kgp_hg19_QC_LD