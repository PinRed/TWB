#!/usr/bin/env bash
set -Eeuo pipefail

############################################################
# F-SNPQC-05
# Within-array batch-effect SNP QC
#
# Input:
#   {ARRAY}_F_SNPQC04_hwe.bed
#   {ARRAY}_F_SNPQC04_hwe.bim
#   {ARRAY}_F_SNPQC04_hwe.fam
#
# Batch mapping:
#   Sample_Batch/TWB2_B0001_samples.keep
#   ...
#   Sample_Batch/TWB2_B0154_samples.keep
#
# Screening rules:
#   Pairwise SNP missingness difference > 0.02
#   Pairwise allele-frequency difference > 0.10
#
# Important:
#   - No additional --geno or --mind filtering.
#   - This script does not perform batch-GWAS.
#   - Threshold flags are candidates for investigation,
#     not proof of technical artifacts.
############################################################

############################################################
# 1. Configuration
############################################################

ARRAY="TWB2"

BFILE="SNP_QC_temp/${ARRAY}_F_SNPQC04_hwe"
BATCH_DIR="../Sample_Batch"
OUT_DIR="${ARRAY}_F_SNPQC05_batch_effect"

FIRST_BATCH=1
LAST_BATCH=154

MISSINGNESS_THRESHOLD=0.02
AF_THRESHOLD=0.10

mkdir -p "${OUT_DIR}"

STATUS="${OUT_DIR}/run_status.txt"
BATCH_LIST="${OUT_DIR}/processed_batches.txt"
SAMPLE_COUNTS="${OUT_DIR}/batch_sample_counts.tsv"

: > "${STATUS}"
: > "${BATCH_LIST}"

echo -e "Batch\tN_samples" > "${SAMPLE_COUNTS}"

############################################################
# 2. Check dependencies and input files
############################################################

for EXT in bed bim fam; do
    if [[ ! -s "${BFILE}.${EXT}" ]]; then
        echo "ERROR: Missing ${BFILE}.${EXT}" >&2
        exit 1
    fi
done

echo "Input: ${BFILE}" | tee -a "${STATUS}"
echo "Array: ${ARRAY}" | tee -a "${STATUS}"
echo "Missingness threshold: ${MISSINGNESS_THRESHOLD}" \
    | tee -a "${STATUS}"
echo "AF threshold: ${AF_THRESHOLD}" | tee -a "${STATUS}"

############################################################
# 3. Record input sample count
############################################################

INPUT_N=$(wc -l < "${BFILE}.fam")

echo "Input samples in F-SNPQC-04: ${INPUT_N}" \
    | tee -a "${STATUS}"

############################################################
# 4. Extract each batch and calculate statistics
############################################################

for ((i=FIRST_BATCH; i<=LAST_BATCH; i++)); do

    printf -v BATCH "%04d" "${i}"

    KEEP="${BATCH_DIR}/${ARRAY}_B${BATCH}_samples.keep"

    if [[ ! -s "${KEEP}" ]]; then
        echo "SKIPPED B${BATCH}: missing/empty keep file" \
            | tee -a "${STATUS}"
        continue
    fi

    BATCH_OUT="${OUT_DIR}/B${BATCH}"
    mkdir -p "${BATCH_OUT}"

    PREFIX="${BATCH_OUT}/${ARRAY}_B${BATCH}"

    echo "=========================================="
    echo "Processing B${BATCH}"
    echo "=========================================="

    ########################################################
    # Step 4.1: Extract frozen F-SNPQC-04 samples by batch
    #
    # Do not apply --geno or --mind here.
    # Preserve variants for batch-effect screening.
    ########################################################

    plink2 \
        --bfile "${BFILE}" \
        --keep "${KEEP}" \
        --make-bed \
        --out "${PREFIX}_base"

    ########################################################
    # Step 4.2: Check output and record sample count
    ########################################################

    if [[ ! -s "${PREFIX}_base.bed" ||
          ! -s "${PREFIX}_base.bim" ||
          ! -s "${PREFIX}_base.fam" ]]; then
        echo "ERROR: Incomplete files for B${BATCH}" >&2
        exit 1
    fi

    N=$(wc -l < "${PREFIX}_base.fam")

    if [[ "${N}" -eq 0 ]]; then
        echo "ERROR: No matched samples in B${BATCH}" >&2
        exit 1
    fi

    echo -e "B${BATCH}\t${N}" >> "${SAMPLE_COUNTS}"

    echo "B${BATCH}: matched samples = ${N}" \
        | tee -a "${STATUS}"

    ########################################################
    # Step 4.3: SNP missingness and allele frequency
    #
    # PLINK2 output:
    #   .vmiss = variant missingness
    #   .afreq = allele frequency
    #
    # No SNP filtering is applied in this step.
    ########################################################

    plink2 \
        --bfile "${PREFIX}_base" \
        --missing \
        --freq \
        --out "${PREFIX}_stats"

    if [[ ! -s "${PREFIX}_stats.vmiss" ||
          ! -s "${PREFIX}_stats.afreq" ]]; then
        echo "ERROR: Missing statistics for B${BATCH}" >&2
        exit 1
    fi

    echo "${BATCH}" >> "${BATCH_LIST}"

    echo "B${BATCH}: statistics completed" \
        | tee -a "${STATUS}"

done

############################################################
# 5. Check batch count
############################################################

N_BATCH=$(wc -l < "${BATCH_LIST}")

echo "Processed batches: ${N_BATCH}" | tee -a "${STATUS}"

if [[ "${N_BATCH}" -lt 2 ]]; then
    echo "ERROR: At least two batches are required." >&2
    exit 1
fi