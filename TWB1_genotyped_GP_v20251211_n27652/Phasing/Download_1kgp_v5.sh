#!/bin/bash
set -Eeuo pipefail

BASE_URL="https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/"
OUTDIR="1kgp_v5"
MAX_ATTEMPTS=10

# Check dependencies
for cmd in wget bcftools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd is not installed" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"
cd "$OUTDIR"

log() {
    echo "[$(date '+%F %T')] $*"
}

get_vcf_name() {
    echo "chr${1}.1kg.phase3.v5a.vcf.gz"
}

validate_vcf() {
    local vcf="$1"

    # Read and decode the entire VCF
    if ! bcftools view -Ov "$vcf" >/dev/null; then
        log "WARNING: VCF full-read check failed: $vcf"
        return 1
    fi

    return 0
}

download_and_validate() {
    local chr="$1"
    local vcf
    local attempt

    vcf="$(get_vcf_name "$chr")"

    # Skip previously validated files
    if [[ -s "$vcf" ]]; then
        log "chr${chr}: existing file found; validating"

        if validate_vcf "$vcf"; then
            log "chr${chr}: existing VCF valid; skipping"
            return 0
        fi

        log "chr${chr}: existing VCF invalid; removing"
        rm -f "$vcf"
    fi

    for ((attempt=1; attempt<=MAX_ATTEMPTS; attempt++)); do

        log "chr${chr}: attempt ${attempt}/${MAX_ATTEMPTS}"

        rm -f "$vcf"

        if ! wget \
            --timeout=60 \
            --tries=3 \
            --waitretry=5 \
            -O "$vcf" \
            "${BASE_URL}${vcf}"; then

            log "WARNING: chr${chr} download failed"
            rm -f "$vcf"
            continue
        fi

        log "chr${chr}: download completed; validating"

        if ! validate_vcf "$vcf"; then
            log "WARNING: chr${chr} validation failed"
            rm -f "$vcf"
            continue
        fi

        log "chr${chr}: download and validation successful"
        return 0
    done

    log "ERROR: chr${chr} failed after ${MAX_ATTEMPTS} attempts"
    return 1
}

for chr in {1..22}; do
    download_and_validate "$chr"
done

log "All chromosomes downloaded and validated successfully"
