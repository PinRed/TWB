#!/bin/bash
set -Eeuo pipefail

ARRAY="TWB1"
INPUT_DIR="../SNP_QC/Harmonization/VCF_Files"
REF_DIR="1kgp_v5"
OUT_DIR="Phased_Impute_VCF"
LOG_DIR="${OUT_DIR}/logs"
BEAGLE="../../Tools/beagle.08Feb22.fa4.jar"

THREADS=4
JAVA_MEM="100g"

mkdir -p "$OUT_DIR" "$LOG_DIR"

FAILED_FILE="${LOG_DIR}/failed_chr.txt"
OOM_FILE="${LOG_DIR}/java_oom_chr.txt"

: > "$FAILED_FILE"
: > "$OOM_FILE"

for chr in {1..22}; do

    echo "=================================="
    echo "Processing chr${chr}"
    echo "Start: $(date)"
    echo "=================================="

    INPUT="${INPUT_DIR}/${ARRAY}_chr${chr}.vcf.gz"
    REF="${REF_DIR}/chr${chr}.EAS.vcf.gz"
    MAP="map/plink.chr${chr}.GRCh37.map"

    PHASE_OUT="${OUT_DIR}/${ARRAY}_chr${chr}_phased"
    IMPUTE_OUT="${OUT_DIR}/${ARRAY}_chr${chr}_imputed"

    PHASE_LOG="${LOG_DIR}/chr${chr}.phase.resource.log"
    IMPUTE_LOG="${LOG_DIR}/chr${chr}.impute.resource.log"

    # ----------------------------
    # Step 1: Phasing
    # ----------------------------

    echo "chr${chr}: PHASING"

    if ! /usr/bin/time -v \
        java -Xms4g -Xmx${JAVA_MEM} \
        -jar "$BEAGLE" \
        gt="$INPUT" \
        ref="$REF" \
        map="$MAP" \
        chrom="$chr" \
        impute=false \
        nthreads="$THREADS" \
        out="$PHASE_OUT" \
        > "$PHASE_LOG" 2>&1
    then
        echo "chr${chr}: PHASING FAILED"
        echo "$chr" >> "$FAILED_FILE"

        if grep -q "OutOfMemoryError" "$PHASE_LOG"; then
            echo "$chr" >> "$OOM_FILE"
        fi

        continue
    fi

    # ----------------------------
    # Step 2: Imputation
    # ----------------------------

    echo "chr${chr}: IMPUTATION"

    if ! /usr/bin/time -v \
        java -Xms4g -Xmx${JAVA_MEM} \
        -jar "$BEAGLE" \
        gt="${PHASE_OUT}.vcf.gz" \
        ref="$REF" \
        map="$MAP" \
        chrom="$chr" \
        impute=true \
        nthreads="$THREADS" \
        out="$IMPUTE_OUT" \
        > "$IMPUTE_LOG" 2>&1
    then
        echo "chr${chr}: IMPUTATION FAILED"
        echo "$chr" >> "$FAILED_FILE"

        if grep -q "OutOfMemoryError" "$IMPUTE_LOG"; then
            echo "$chr" >> "$OOM_FILE"
        fi

        continue
    fi

    echo "chr${chr}: SUCCESS"
    echo "End: $(date)"
    echo

done

echo "=================================="
echo "FINAL SUMMARY"
echo "=================================="

echo "Failed chromosomes:"
if [[ -s "$FAILED_FILE" ]]; then
    cat "$FAILED_FILE"
else
    echo "None"
fi

echo
echo "Java OOM chromosomes:"
if [[ -s "$OOM_FILE" ]]; then
    cat "$OOM_FILE"
else
    echo "None"
fi

if [[ -s "$FAILED_FILE" ]]; then
    exit 1
fi

echo "All chromosomes completed successfully."