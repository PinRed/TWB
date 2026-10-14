#!/usr/bin/env bash
# ========================================
# Script: 4-1_hg19_to_hg38.sh
# Purpose: Convert PLINK dataset from hg19 to hg38 using LiftOver
# ========================================

set -euo pipefail

# hg19 轉 hg38
# Tools : liftover

########################################
# Settings
########################################

INPUT="../King/TWB1_S05"
OUTPUT="TWB1_S05_hg38"

CHAIN="../../hg19ToHg38.over.chain.gz"
LIFTOVER="liftOver"

# PLINK command
PLINK="plink"


########################################
# Check files
########################################

echo "========================================"
echo "PLINK hg19 -> hg38 Liftover"
echo "========================================"

if [[ ! -f "${INPUT}.bed" ]]; then
    echo "ERROR: ${INPUT}.bed not found"
    exit 1
fi

if [[ ! -f "${INPUT}.bim" ]]; then
    echo "ERROR: ${INPUT}.bim not found"
    exit 1
fi

if [[ ! -f "${INPUT}.fam" ]]; then
    echo "ERROR: ${INPUT}.fam not found"
    exit 1
fi

if [[ ! -f "${CHAIN}" ]]; then
    echo "ERROR: ${CHAIN} not found"
    echo ""
    echo "Download with:"
    echo "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
    exit 1
fi

if ! command -v "${LIFTOVER}" >/dev/null 2>&1; then
    echo "ERROR: liftOver command not found"
    exit 1
fi

if ! command -v "${PLINK}" >/dev/null 2>&1; then
    echo "ERROR: plink command not found"
    exit 1
fi


########################################
# Temporary files
########################################

TMP="liftover_tmp"

rm -rf "${TMP}"
mkdir -p "${TMP}"

BIM="${INPUT}.bim"

INPUT_BED="${TMP}/input_hg19.bed"
OUTPUT_BED="${TMP}/output_hg38.bed"
UNMAPPED_BED="${TMP}/unmapped.bed"

UPDATE_MAP="${TMP}/update_map.txt"
KEEP_SNPS="${TMP}/mapped_snps.txt"
EXCLUDE_SNPS="${TMP}/unmapped_snps.txt"


########################################
# Step 1
# Check duplicate SNP IDs
########################################

echo ""
echo "[1/7] Checking SNP IDs..."

DUPLICATES=$(awk '{print $2}' "${BIM}" | sort | uniq -d | head)

if [[ -n "${DUPLICATES}" ]]; then
    echo ""
    echo "ERROR: Duplicate SNP IDs detected:"
    echo "${DUPLICATES}"
    echo ""
    echo "Please resolve duplicate SNP IDs before running this script."
    exit 1
fi

echo "No duplicate SNP IDs found."


########################################
# Step 2
# Convert BIM -> BED
########################################

echo ""
echo "[2/7] Converting BIM to BED..."

awk 'BEGIN{OFS="\t"}
{
    chr=$1

    # LiftOver expects chromosome names such as chr1
    if (chr !~ /^chr/) {
        chr="chr" chr
    }

    # BED is 0-based, half-open
    start=$4-1
    end=$4

    # Use SNP ID as BED name
    print chr, start, end, $2
}' "${BIM}" > "${INPUT_BED}"

echo "Created:"
echo "  ${INPUT_BED}"


########################################
# Step 3
# Run LiftOver
########################################

echo ""
echo "[3/7] Running hg19 -> hg38 LiftOver..."

"${LIFTOVER}" \
    "${INPUT_BED}" \
    "${CHAIN}" \
    "${OUTPUT_BED}" \
    "${UNMAPPED_BED}"

echo ""
echo "LiftOver finished."


########################################
# Step 4
# Create mapped / unmapped SNP lists
########################################

echo ""
echo "[4/7] Creating SNP lists..."

# Successfully mapped SNPs
awk 'BEGIN{OFS="\t"} {print $4}' \
    "${OUTPUT_BED}" > "${KEEP_SNPS}"

# Unmapped SNPs
awk '
    /^#/ {next}
    NF >= 4 {print $4}
' "${UNMAPPED_BED}" > "${EXCLUDE_SNPS}"


MAPPED=$(wc -l < "${KEEP_SNPS}")
UNMAPPED=$(wc -l < "${EXCLUDE_SNPS}")
TOTAL=$(wc -l < "${BIM}")

echo "Total SNPs   : ${TOTAL}"
echo "Mapped SNPs  : ${MAPPED}"
echo "Unmapped SNPs: ${UNMAPPED}"


########################################
# Step 5
# Create PLINK update-map file
########################################

echo ""
echo "[5/7] Creating hg38 position file..."

# update_map format:
# SNP_ID  CHR  BP

awk 'BEGIN{OFS="\t"}
{
    print $4, $3
}' "${OUTPUT_BED}" > "${UPDATE_MAP}"

echo "Created:"
echo "  ${UPDATE_MAP}"


########################################
# Step 6
# Generate new PLINK dataset
########################################

echo ""
echo "[6/7] Creating new PLINK dataset..."

# First remove SNPs that could not be mapped.
# Then update chromosome and position.

plink \
    --bfile "${INPUT}" \
    --exclude "${EXCLUDE_SNPS}" \
    --update-map "${UPDATE_MAP}" \
    --make-bed \
    --out "${TMP}/${OUTPUT}"

plink2 \
    --bfile "${TMP}/${OUTPUT}" \
    --set-all-var-ids '@:#:$r:$a' \
    --make-bed \
    --out "${OUTPUT}"

########################################
# Step 7
# Final check
########################################

echo ""
echo "[7/7] Checking output..."

if [[ ! -f "${OUTPUT}.bed" ]]; then
    echo "ERROR: ${OUTPUT}.bed was not created."
    exit 1
fi

if [[ ! -f "${OUTPUT}.bim" ]]; then
    echo "ERROR: ${OUTPUT}.bim was not created."
    exit 1
fi

if [[ ! -f "${OUTPUT}.fam" ]]; then
    echo "ERROR: ${OUTPUT}.fam was not created."
    exit 1
fi

echo ""
echo "========================================"
echo "DONE"
echo "========================================"

echo ""
echo "Output files:"
echo "  ${OUTPUT}.bed"
echo "  ${OUTPUT}.bim"
echo "  ${OUTPUT}.fam"

echo ""
echo "SNP count:"
wc -l "${OUTPUT}.bim"

echo ""
echo "First 5 variants:"
head -5 "${OUTPUT}.bim"

echo ""
echo "Temporary files:"
echo "  ${TMP}/"

echo ""
echo "========================================"


rm -rf "${TMP}"