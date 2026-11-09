
#!/bin/bash
set -Eeuo pipefail

# ==========================================
# TWB1 Eagle2 + Minimac4
# GRCh37 / 1000G EAS
# ==========================================

ARRAY="TWB1"

INPUT_DIR="../SNP_QC/Harmonization/VCF_Files"
REF_DIR="1kgp_v5/EAS_filtered"
REF_DIR_MSAV="1kgp_v5/Minimac_ReHap"
MAP_DIR="/SATA2TB/Taiwan_Biobank/Tools/Eagle_v2.4.1/tables"

OUT_DIR="Phased_Impute_VCF"
LOG_DIR="${OUT_DIR}/logs"

EAGLE="../../Tools/Eagle_v2.4.1/eagle"
MINIMAC4="Minimac4"

THREADS=16

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

    INPUT="${INPUT_DIR}/${ARRAY}_chr${chr}.vcf.gz"

    REF_VCF="${REF_DIR}/chr${chr}.EAS.vcf.gz"
    REF_MSAV="${REF_DIR_MSAV}/chr${chr}.EAS.msav"

    # Eagle-compatible GRCh37 genetic map
    MAP="${MAP_DIR}/genetic_map_hg19_withX.txt.gz"

    PHASE_OUT="${OUT_DIR}/${ARRAY}_chr${chr}_phased"
    IMPUTE_OUT="${OUT_DIR}/${ARRAY}_chr${chr}_imputed"

    PHASE_VCF="${PHASE_OUT}.vcf.gz"

    PHASE_LOG="${LOG_DIR}/chr${chr}.phase.resource.log"
    IMPUTE_LOG="${LOG_DIR}/chr${chr}.impute.resource.log"

    # --------------------------------------
    # Input validation
    # --------------------------------------

    if [[ ! -s "$INPUT" ||
          ! -s "$REF_VCF" ||
          ! -s "$REF_MSAV" ||
          ! -s "$MAP" ]]; then

        echo "chr${chr}: Missing input file"
        echo "$chr" >> "$FAILED_FILE"
        continue
    fi

    # --------------------------------------
    # Step 1: Eagle Phasing
    # --------------------------------------

    echo "chr${chr}: EAGLE PHASING"

    if ! /usr/bin/time -v \
        "$EAGLE" \
        --vcfRef "$REF_VCF" \
        --vcfTarget "$INPUT" \
        --geneticMapFile "$MAP" \
        --chrom "$chr" \
        --allowRefAltSwap \
        --vcfOutFormat z \
        --numThreads "$THREADS" \
        --outPrefix "$PHASE_OUT" \
        > "$PHASE_LOG" 2>&1
    then

        echo "chr${chr}: PHASING FAILED"
        echo "$chr" >> "$FAILED_FILE"

        if grep -Eqi \
            "out of memory|bad_alloc|cannot allocate memory" \
            "$PHASE_LOG"; then
            echo "$chr:phasing" >> "$OOM_FILE"
        fi

        continue
    fi

    # --------------------------------------
    # Step 2: Minimac4 Imputation
    # --------------------------------------

    echo "chr${chr}: MINIMAC4 IMPUTATION"

    if ! /usr/bin/time -v \
        "$MINIMAC4" \
        --refHaps "$REF_MSAV=" \
        --haps "$PHASE_VCF" \
        --window 500000 \
        --prefix "$IMPUTE_OUT" \
        --cpus "$THREADS" \
        --chr "$chr" \
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
