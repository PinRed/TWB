#!/bin/bash
set -Eeuo pipefail

BASE_URL="https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/"
OUTDIR="1kgp_v5"
MAX_ATTEMPTS=10

mkdir -p "$OUTDIR"
cd "$OUTDIR"

log() {
    echo "[$(date '+%F %T')] $*"
}

get_vcf_name() {
    echo "chr${1}.1kg.phase3.v5a.vcf.gz"
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

        return 0
    done

    log "ERROR: chr${chr} failed after ${MAX_ATTEMPTS} attempts"
    return 1
}

for chr in {1..22}; do
    download_and_validate "$chr"
done
