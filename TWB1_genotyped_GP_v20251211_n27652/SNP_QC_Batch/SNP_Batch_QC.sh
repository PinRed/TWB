
#!/bin/bash

############################################################
# TWB1 Batch SNP QC
# PLINK2 + PLINK binary fileset (.bed/.bim/.fam)
#
# Input:
#   TWB1_S03_QC.bed
#   TWB1_S03_QC.bim
#   TWB1_S03_QC.fam
#
# Sample lists:
#   Sample_Batch/TWB1_B0001_samples.keep
#   ...
#   Sample_Batch/TWB1_B0154_samples.keep
#
# QC:
#   1. SNP missing rate > 10%       --geno 0.10
#   2. Individual missing rate >10% --mind 0.10
#   3. SNP missing rate > 5%        --geno 0.05
#   4. Calculate missingness         --missing
#   5. Calculate allele frequency    --freq
#   6. SNP intersection across all batches
############################################################


############################################################
# Settings
############################################################

BFILE="../Primary_Sample_QC/TWB1_S03_QC"
BATCH_DIR="../Sample_Batch"
OUT_DIR="Batch_SNP_QC"

mkdir -p "${OUT_DIR}"


############################################################
# Check input files
############################################################

for EXT in bed bim fam
do
    if [[ ! -s "${BFILE}.${EXT}" ]]; then
        echo "ERROR: Missing ${BFILE}.${EXT}"
        exit 1
    fi
done


############################################################
# 1. Batch-level SNP QC
############################################################

for BATCH in {1..14}
do

    KEEP="${BATCH_DIR}/TWB1_B${BATCH}_samples.keep"

    if [[ ! -s "${KEEP}" ]]; then
        echo "Skipping batch ${BATCH}: ${KEEP} missing or empty"
        continue
    fi

    echo "=========================================="
    echo "Processing Batch ${BATCH}"
    echo "=========================================="

    BATCH_OUT="${OUT_DIR}/B${BATCH}"

    mkdir -p "${BATCH_OUT}"


    ########################################################
    # Step 1
    # Keep samples belonging to this batch
    # Remove SNPs with missing rate > 10%
    ########################################################

    plink2 \
        --bfile "${BFILE}" \
        --keep "${KEEP}" \
        --geno 0.10 \
        --make-bed \
        --out "${BATCH_OUT}/TWB1_B${BATCH}_geno01"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Batch ${BATCH} failed at Step 1"
        exit 1
    fi


    ########################################################
    # Step 2
    # Remove individuals with missing rate > 10%
    ########################################################

    plink2 \
        --bfile "${BATCH_OUT}/TWB1_B${BATCH}_geno01" \
        --mind 0.10 \
        --make-bed \
        --out "${BATCH_OUT}/TWB1_B${BATCH}_geno01_mind01"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Batch ${BATCH} failed at Step 2"
        exit 1
    fi


    ########################################################
    # Step 3
    # Remove SNPs with missing rate > 5%
    # Calculate SNP missingness and allele frequency
    ########################################################

    plink2 \
        --bfile "${BATCH_OUT}/TWB1_B${BATCH}_geno01_mind01" \
        --geno 0.05 \
        --make-bed \
        --missing \
        --freq \
        --out "${BATCH_OUT}/TWB1_B${BATCH}_geno01_mind01_geno005"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Batch ${BATCH} failed at Step 3"
        exit 1
    fi


    ########################################################
    # Check final output
    ########################################################

    if [[ ! -s "${BATCH_OUT}/TWB1_B${BATCH}_geno01_mind01_geno005.bim" ]]; then
        echo "ERROR: Missing final BIM for Batch ${BATCH}"
        exit 1
    fi

    echo "Batch ${BATCH} completed."

done


############################################################
# 2. SNP intersection across all batches
############################################################

echo ""
echo "=========================================="
echo "Calculating SNP intersection"
echo "=========================================="


############################################################
# Find all final BIM files
############################################################

mapfile -t BIM_LIST < <(
    find "${OUT_DIR}" \
        -mindepth 2 \
        -maxdepth 2 \
        -name '*_geno01_mind01_geno005.bim' \
        | sort
)


if [[ ${#BIM_LIST[@]} -eq 0 ]]; then
    echo "ERROR: No final BIM files found."
    exit 1
fi


echo "Number of batches found: ${#BIM_LIST[@]}"


############################################################
# 3. Start with first batch
############################################################

FIRST_BIM="${BIM_LIST[0]}"

echo "Starting BIM:"
echo "${FIRST_BIM}"

sort "${FIRST_BIM}" > a


############################################################
# 4. Intersection of SNPs across all batches
############################################################

: > snp_intersect.log

for BIM in "${BIM_LIST[@]}"
do

    echo "Intersecting:"
    echo "${BIM}"

    sort "${BIM}" > batch_bim

    comm -12 a batch_bim > b

    mv b a

    COUNT=$(wc -l < a)

    echo "${BIM} ${COUNT}" >> snp_intersect.log

done


rm -f batch_bim


############################################################
# 5. Final intersection SNP list
############################################################

mv a TWB1_SNP_Batch_QC.bim

cut -f2 TWB1_SNP_Batch_QC.bim > TWB1_SNP_Batch_QC.id


############################################################
# Summary
############################################################

FINAL_COUNT=$(wc -l < TWB1_SNP_Batch_QC.id)

echo ""
echo "=========================================="
echo "SNP intersection completed"
echo "Final SNP count: ${FINAL_COUNT}"
echo "=========================================="

echo "Output:"
echo "  TWB1_SNP_Batch_QC.bim"
echo "  TWB1_SNP_Batch_QC.id"
echo "  snp_intersect.log"

rm -rf "${OUT_DIR}" 


####################################################
# 抓共同 SNP ID
# 1. 抓出 SNP ID 並排序
awk '{print $2}' SNP_QC/TWB1_S03_all_QC_clean.bim | sort > list_a.txt
awk '{print $2}' TWB1_SNP_Batch_QC.bim | sort > list_b.txt

# 2. 取兩者的交集
comm -12 list_a.txt list_b.txt > intersect_snps.id

# 3. 清除暫存檔
rm -f list_a.txt list_b.txt