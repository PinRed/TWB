
#!/bin/bash
set -Eeuo pipefail

# ==========================================
# TWB2 Minimac4 Imputation
# GRCh38 / 1000G EAS
# ==========================================

ARRAY="TWB2"

REF_DIR_MSAV="1kgp_v5/Minimac_ReHap"

OUT_DIR="Phased_Impute_VCF"
LOG_DIR="${OUT_DIR}/logs"

MINIMAC4="minimac4"
THREADS=24

mkdir -p "$OUT_DIR" "$LOG_DIR"

FAILED_FILE="${LOG_DIR}/failed_chr.txt"
OOM_FILE="${LOG_DIR}/oom_chr.txt"
SUCCESS_FILE="${LOG_DIR}/success_chr.txt"

: > "$FAILED_FILE"
: > "$OOM_FILE"
: > "$SUCCESS_FILE"

for chr in {1..22}; do

    echo "=================================="
    echo "Processing chr${chr}"
    echo "Start: $(date)"
    echo "=================================="

    REF_MSAV="${REF_DIR_MSAV}/chr${chr}.EAS.msav"

    PHASE_VCF="${OUT_DIR}/${ARRAY}_chr${chr}_phased.vcf.gz"

    IMPUTE_OUT="${OUT_DIR}/${ARRAY}_chr${chr}_imputed"
    IMPUTE_VCF="${IMPUTE_OUT}.dose.vcf.gz"

    IMPUTE_LOG="${LOG_DIR}/chr${chr}.impute.resource.log"

    # --------------------------------------
    # Check existing imputed VCF
    # --------------------------------------

    if [[ -s "$IMPUTE_VCF" ]] &&
       gzip -t "$IMPUTE_VCF" 2>/dev/null; then

        echo "chr${chr}: IMPUTED VCF EXISTS, SKIPPING MINIMAC4"
        echo "$chr" >> "$SUCCESS_FILE"
        continue
    fi

    # --------------------------------------
    # Input validation
    # --------------------------------------

    if [[ ! -s "$PHASE_VCF" ||
          ! -s "$REF_MSAV" ]]; then

        echo "chr${chr}: Missing input file"
        echo "$chr" >> "$FAILED_FILE"
        continue
    fi

    # --------------------------------------
    # Minimac4 Imputation
    # --------------------------------------

    echo "chr${chr}: MINIMAC4 IMPUTATION"

    if ! /usr/bin/time -v \
        "$MINIMAC4" \
        --refHaps "$REF_MSAV" \
        --haps "$PHASE_VCF" \
        --window 500000 \
        --prefix "$IMPUTE_OUT" \
        --cpus "$THREADS" \
        --chr "chr${chr}" \
        --noPhoneHome \
        --format GT,DS,GP \
        --allTypedSites \
        --meta \
        --minRatio 0.00001 \
        > "$IMPUTE_LOG" 2>&1
    then

        echo "chr${chr}: IMPUTATION FAILED"
        echo "$chr" >> "$FAILED_FILE"

        if grep -Eqi \
            "out of memory|bad_alloc|cannot allocate memory" \
            "$IMPUTE_LOG"; then
            echo "$chr:imputation" >> "$OOM_FILE"
        fi

        continue
    fi

    # --------------------------------------
    # Validate imputed output
    # --------------------------------------

    if [[ ! -s "$IMPUTE_VCF" ]] ||
       ! gzip -t "$IMPUTE_VCF" 2>/dev/null; then

        echo "chr${chr}: IMPUTED VCF INVALID"
        echo "$chr" >> "$FAILED_FILE"
        continue
    fi

    echo "$chr" >> "$SUCCESS_FILE"

    echo "chr${chr}: SUCCESS"
    echo "End: $(date)"
    echo

done

echo "=================================="
echo "FINAL SUMMARY"
echo "=================================="

echo "Successful chromosomes:"
if [[ -s "$SUCCESS_FILE" ]]; then
    cat "$SUCCESS_FILE"
else
    echo "None"
fi

echo
echo "Failed chromosomes:"
if [[ -s "$FAILED_FILE" ]]; then
    sort -nu "$FAILED_FILE"
else
    echo "None"
fi

echo
echo "Memory errors:"
if [[ -s "$OOM_FILE" ]]; then
    cat "$OOM_FILE"
else
    echo "None"
fi

if [[ -s "$FAILED_FILE" ]]; then
    exit 1
fi

echo "All chromosomes completed successfully."
