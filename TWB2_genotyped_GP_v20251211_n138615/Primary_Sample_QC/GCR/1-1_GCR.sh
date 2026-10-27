#!/bin/bash
# ==============================================================================
# Pipeline Step: TWB2 Batch Sample Missing Rate Calculation
# Target Dataset: TWB2 cohort 
# Inputs:
#   - TWB2.hg38.GP.20251211.{bed,bim,fam} (TWB2 genotyped dataset)
#   - TWB2_B{0001..0154}_samples.keep (List of samples for each batch)
# Outputs:
#   - Missing/TWB2_Batch{0001..0154}_missing.{imiss,imiss.log} (Missing rate for each batch)
#   - TWB2_S01_GCR.{bed,bim,fam} (TWB2 dataset after removing samples failing GCR QC)
# Do:
#   1. Calculate missing rate for each batch of TWB2 samples
#   2. Remove samples failing GCR QC and create a new dataset
# ==============================================================================


PREFIX="../../TWB2.hg38.GP.20251211"
# TWB 2 Batch sample missing rate 
for BATCH in {0001..0154}; do
    KEEP="../../Sample_Batch/TWB2_B${BATCH}_samples.keep"

    if [[ ! -s "$KEEP" ]]; then
        echo "Skipping batch ${BATCH}: ${KEEP} missing or empty"
        continue
    fi

    plink2 \
        --bfile "$PREFIX" \
        --keep "$KEEP" \
        --missing \
        --out "Missing/TWB2_Batch${BATCH}_missing"
done

# TWB_Sample_QC.R
# output TWB2_S01_GCR_exclude.txt

plink2 --bfile "$PREFIX" --remove TWB2_S01_GCR_exclude.txt --make-bed --out TWB2_S01_GCR